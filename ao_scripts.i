
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

require,yoga_ao_top+"/yoga_ao.i";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

YOGA_AO_PARPATH = YOGA_AO_SAVEPATH+"/par/";
mkdirp,YOGA_AO_PARPATH;

func script_atmos(filename)
{
  extern y_geom,y_tel,y_loop,y_atmos;
  extern g_atmos,g_target;
  extern ipupil;

  if (filename == []) filename = YOGA_AO_PARPATH+"atmos_1layer.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"atmos_12layers.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  read_parfile,filename;

  if (y_loop.niter == []) y_loop.niter = 1000;
  
  //here we just simulate atmosphere so pupdiam is fixed
  // arbitrarily.
  pupdiam = long(32*y_tel.diam/4.0);

  write,"Doing atmos inits on the GPU";

  geom_init,pupdiam;

  atmos_init;
  
  write,"... Done !";
  
  write,"Creating targets on the GPU";

  target_init;

  write,"... Done !";

  write,format="Starting loop on : %d iterations\n",y_loop.niter;
  
  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    //tinter = tic();
    move_sky,g_atmos,g_target;
    if (cc % 10 == 0) {
      time_move = tac(mytime)/cc;
      write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
    }
  }
  write,"... Done !";
  write,format="Average atmos move gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  
  write,"Doing atmos inits on the CPU";

  pupixsize = y_tel.diam / y_geom.pupdiam;
  
  cpu_r0 = y_atmos.r0;

  screen_cpu = create_screen(cpu_r0,pupixsize,y_geom.pupdiam,A,B,ist);

  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    //tinter = tic();
    screen_cpu = extrude(screen_cpu,float(cpu_r0 /pupixsize), A, B, ist);
    if (cc % 10 == 0) {
      time_move = tac(mytime)/cc;
      write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
    }
  }
  write,"... Done !";
  write,format="Average extrude cpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  
  // now that we have both gpu and cpu inits, compare screen
  // need to compute phase structure function
  tmp = array(float,[2,2*y_geom.pupdiam,2*y_geom.pupdiam]);
  tmp(1:y_geom.pupdiam,1:y_geom.pupdiam)  = *y_geom._spupil;
  den  = (fft(fft(tmp,1)*conj(fft(tmp,1)),-1)).re;
  dphi_cpu = dphi_gpu1 = dphi_gpu2 = 0.0f;
  
  write,"Comparing CPU vs GPU screens";

  niter = 100000;
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(cpu_r0 /pupixsize), A, B, ist);
     dphi_cpu   += calc_dphi(screen_cpu(1:y_geom.pupdiam,1:y_geom.pupdiam),*y_geom._spupil,den);
     move_sky,g_atmos,g_target;
     //screen_gpu1 = get_tscreen(g_atmos,0.)(1:y_geom.pupdiam,1:y_geom.pupdiam);
     target_atmostrace,g_target,0,g_atmos;
     screen_gpu2 = target_getphase(g_target,0);
     //dphi_gpu1   += calc_dphif(screen_gpu1,*y_geom._spupil,den,conjftpup);
     dphi_gpu2   += calc_dphi(screen_gpu2,*y_geom._spupil,den);
     if (cc % 10 == 0) {
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
     }
  }
  
  write,"... Done !";

  fma;
  plg,6.88*(indgen(y_geom.pupdiam)/float(cpu_r0 /pupixsize))^(5./3.),marks=0,width=4;
  plg,circavg_quad(dphi_cpu(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter,color="red",marks=0,width=4;
  plg,circavg_quad(dphi_gpu2(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter,color="green",marks=0,width=4;
  pltitle,"<(phase-phase(1,1))^2^ > , GPU (g) vs CPU (r)"
}

func script_atmoscpu(void)
{
  niter = 100000;
  pupdiam = 64;
  
  spupil = make_pupil(pupdiam,pupdiam,xc=pupdiam/2+0.5,yc=pupdiam/2+0.5,cobs=0.1)
  tmp = array(float,[2,2*pupdiam,2*pupdiam]);
  tmp(1:pupdiam,1:pupdiam)  = spupil;
  conjftpup = conj(fft(tmp,1));
  den  = (fft(fft(tmp,1)*conj(fft(tmp,1)),-1)).re;
  dphi_cpu0 = dphi_cpu1 = dphi_cpu2 = dphi_cpu3 = 0.0f;

  screen_cpu = create_screen(0.16,0.0313725,pupdiam,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    //dphi_cpu0   += calc_dphif(screen_cpu,spupil,den,conjftpup);
      dphi_cpu0   += calc_dphi(screen_cpu,spupil,den);
    time_move = tac(mytime)/cc;
    write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam+16,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    //dphi_cpu1   += calc_dphif(screen_cpu(1:pupdiam,1:pupdiam),spupil,den,conjftpup);
    dphi_cpu1   += calc_dphi(screen_cpu(1:pupdiam,1:pupdiam),spupil,den);
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam*2,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    //dphi_cpu2   += calc_dphif(screen_cpu(1:pupdiam,1:pupdiam),spupil,den,conjftpup);
    dphi_cpu2   += calc_dphi(screen_cpu(1:pupdiam,1:pupdiam),spupil,den);
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }

  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam*4,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
      dphi_cpu3   += calc_dphi(screen_cpu(1:pupdiam,1:pupdiam),spupil,den);
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  plg,6.88*(indgen(pupdiam)/5.1)^(5./3.),marks=0,width=4;
  plg,circavg_quad(dphi_cpu0(1:pupdiam,1:pupdiam))/niter,color="red",marks=0,width=4;
  plg,circavg_quad(dphi_cpu1(1:pupdiam,1:pupdiam))/niter,color="blue",marks=0,width=4;
  plg,circavg_quad(dphi_cpu2(1:pupdiam,1:pupdiam))/niter,color="green",marks=0,width=4;
  plg,circavg_quad(dphi_cpu3(1:pupdiam,1:pupdiam))/niter,color="yellow",marks=0,width=4;
  error;
}
