
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
  dphi_cpu = dphi_gpu1 = dphi_gpu2 = 0.0f;
  
  write,"Comparing CPU vs GPU screens";

  niter = 100000;
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(cpu_r0 /pupixsize), A, B, ist);
     dphi_cpu   += calc_dphis(screen_cpu(1:y_geom.pupdiam,1:y_geom.pupdiam));
     move_sky,g_atmos,g_target;
     target_atmostrace,g_target,0,g_atmos;
     screen_gpu2 = target_getphase(g_target,0);
     dphi_gpu2   += calc_dphis(screen_gpu2);
     if (cc % 10 == 0) {
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
     }
  }
  
  write,"... Done !";

  fma;
  circa = circavg_quad(dphi_gpu2(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter;
  plg,6.88*(indgen(numberof(circa))/float(cpu_r0 /pupixsize))^(5./3.),marks=0,width=4;
  plg,circavg_quad(dphi_cpu(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter,color="red",marks=0,width=4;
  plg,circa,color="green",marks=0,width=4;
  pltitle,"<(phase-phase(1,1))^2^ > , GPU (g) vs CPU (r)"
}

func script_atmosgpu(filename)
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
  
  
  // now that we have both gpu and cpu inits, compare screen
  // need to compute phase structure function
  dphi_gpu = 0.0f;
  
  write,"Comparing CPU vs GPU screens";

  niter = 100000;
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
     move_sky,g_atmos,g_target;
     target_atmostrace,g_target,0,g_atmos;
     screen_gpu = target_getphase(g_target,0);
     dphi_gpu   += calc_dphis(screen_gpu);
     if (cc % 10 == 0) {
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
     }
  }
  
  write,"... Done !";

  fma;
  pupixsize = y_tel.diam / y_geom.pupdiam;
  circa = circavg_quad(dphi_gpu(1:y_geom.pupdiam,1:y_geom.pupdiam))/niter;
  plg,6.88*(indgen(numberof(circa))/float(y_atmos.r0 /pupixsize))^(5./3.),marks=0,width=4;
  plg,circa,color="green",marks=0,width=4;
  pltitle,"<(phase-phase(1,1))^2^ > , GPU (g) vs CPU (r)"
}

func script_atmoscpu(void)
{
  niter = 100000;
  pupdiam = 64;
  
  dphi_cpu0 = dphi_cpu1 = dphi_cpu2 = dphi_cpu3 = 0.0f;

  screen_cpu = create_screen(0.16,0.0313725,pupdiam,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
      dphi_cpu0   += calc_dphis(screen_cpu);
    time_move = tac(mytime)/cc;
    write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam+16,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    dphi_cpu1   += calc_dphis(screen_cpu(1:pupdiam,1:pupdiam));
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam*2,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
    dphi_cpu2   += calc_dphis(screen_cpu(1:pupdiam,1:pupdiam));
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }

  
  A = B = ist = 0.;
  
  screen_cpu = create_screen(0.16,0.0313725,pupdiam*4,A,B,ist);
  mytime = tic();
  for (cc=1;cc<=niter;cc++) {
    screen_cpu = extrude(screen_cpu,float(0.16 /0.0313725), A, B, ist);
      dphi_cpu3   += calc_dphis(screen_cpu(1:pupdiam,1:pupdiam));
       time_move = tac(mytime)/cc;
       write,format="\r Estimated remaining time : %.2f s",(niter - cc)*time_move;
  }
  
  plg,6.88*(indgen(pupdiam)/5.1)^(5./3.),marks=0,width=4;
  plg,circavg_quad(dphi_cpu0(1:pupdiam,1:pupdiam))/niter,color="red",marks=0,width=4;
  plg,circavg_quad(dphi_cpu1(1:pupdiam,1:pupdiam))/niter,color="blue",marks=0,width=4;
  plg,circavg_quad(dphi_cpu2(1:pupdiam,1:pupdiam))/niter,color="green",marks=0,width=4;
  plg,circavg_quad(dphi_cpu3(1:pupdiam,1:pupdiam))/niter,color="yellow",marks=0,width=4;
}

func script_wfs(filename,typeslp)
{

  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (filename == []) filename = YOGA_AO_PARPATH+"atmos_1layer.par";
  //if (filename == []) filename = YOGA_AO_PARPATH+"atmos_12layers.par";

  if ((!(fileExist(filename))) && (!(fileExist(YOGA_AO_PARPATH+filename))))
    error,"could not find"+filename;
  
  if (!(fileExist(filename)))
    filename = YOGA_AO_PARPATH+filename;

  read_parfile,filename;

  if (typeslp == []) typeslp=0;
  
  if (y_loop.niter == []) y_loop.niter = 1000;

  wfs_init;

  atmos_init;
 
  target_init;

  write,"... Done with inits !";

  write,format="Starting loop on : %d iterations\n",y_loop.niter;
  
  mytime = tic();
  for (cc=1;cc<=y_loop.niter;cc++) {
    //tinter = tic();
    move_sky,g_atmos,g_target;
    if ((y_wfs != []) && (g_wfs != [])) {
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,g_atmos;
        sensors_compimg,g_wfs,i-1;
        if (typeslp > 0) {
          if (typeslp == 1) slopes_geom,g_wfs,i-1,0;
          if (typeslp == 2) slopes_geom,g_wfs,i-1,1;
          if (typeslp == 3) sensors_compslopes,g_wfs,i-1;
          if (typeslp == 4) sensors_compslopes,g_wfs,i-1,1,10;
          if (typeslp == 5) sensors_compslopes,g_wfs,i-1,2,100.;
        }
      }
    }
    if (cc % 10 == 0) {
      time_move = tac(mytime)/cc;
      write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
    }
  }
  write,"... Done !";
  write,format="Average wfs gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  
}


func check_centroiding(void)
{
  // check geom slopes
  slopes_geom,g_wfs,0,0;
  slp=sensors_getslopes(g_wfs,0);
  res=reform(sensors_getdata(g_wfs,0,"phase")(*y_wfs(1)._phasemap+1),y_wfs(1)._pdiam,y_wfs(1)._pdiam,y_wfs(1)._nvalid)-(*y_wfs(1)._halfxy)(,,-:1:y_wfs(1)._nvalid);
  slpx = res(0,avg,)-res(1,avg,);
  slpy = res(avg,0,)-res(avg,1,);
  geom=_(slpx,slpy);
  geom/slp;

  slopes_geom,g_wfs,0,0;
  slp=sensors_getslopes(g_wfs,0);
  slopes_geom,g_wfs,0,1;
  slp1=sensors_getslopes(g_wfs,0);slp/slp1;

  // check centroids
  sensors_compslopes,g_wfs,0;slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  slpx=(res2*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  slpy=(res2*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  centro=_(slpx,slpy);
  centro/slp2;

  // check sortmax
  nmax = 16;
  sensors_initnmax,g_wfs,0,nmax;
  sensors_getnmax,g_wfs,0;
  res1=sensors_getdata(g_wfs,0,"validp");
  res1(,0);
  res=sensors_getdata(g_wfs,0,"bincube");
  res(,,0)(*)(sort(res(*,0)))(::-1)(1:nmax);
  res1=sensors_getdata(g_wfs,0,"validi");
  res1(,0);
  (sort(res(*,0)))(::-1)(1:nmax)-1;

  // check centroids nmax
  nmax=10;
  sensors_compslopes,g_wfs,0,1,nmax;
  slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  for (i=1;i<=dimsof(res2)(4);i++) res2(*,i)((sort(res2(*,i)))(::-1)(nmax+1:))=0.;
  slpx=(res2*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  slpy=(res2*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  centro=_(slpx,slpy);
  centro/slp2;

  // check thresholded centroids 
  thresh=2000;
  sensors_compslopes,g_wfs,0,2,thresh;
  slp2=sensors_getslopes(g_wfs,0);
  res2=sensors_getdata(g_wfs,0,"bincube");
  tmp=indices(y_wfs(1).npix);
  for (i=1;i<=dimsof(res2)(4);i++) {
    nn = where(res2(*,i) < thresh);
    if (numberof(nn) > 0) res2(*,i)(nn) *= 0.;
  }
  slpx=(res2*tmp(,,1)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  slpy=(res2*tmp(,,2)(,,-:1:y_wfs(1)._nvalid))(*,)(sum,)/res2(*,)(sum,);
  centro=_(slpx,slpy);
  centro/slp2;

}
/*

tesla 2070
  canaray_elt : 
Average wfs gpu time : 0.1088 s => 9.19118 Hz
  canaray_elt + centroid : 
Average wfs gpu time : 0.1340 s => 7.46269 Hz
  1wfs80x80_1layer :
Average wfs gpu time : 0.0253 s => 39.5257 Hz

gtx 480
  1wfs80x80_1layer :
Average wfs gpu time : 0.0154 s => 64.9351 Hz
  1wfs80x80_1layer + centroid :
Average wfs gpu time : 0.0200 s => 50 Hz
  canary_ngs :
Average wfs gpu time : 0.0016 s => 625 Hz
  canary_ngs + centroid :
Average wfs gpu time : 0.0018 s => 555.556 Hz
  1wfs16x16_1layer :
Average wfs gpu time : 0.0007 s => 1428.57 Hz
  1wfs16x16_1layer + centroid :
Average wfs gpu time : 0.0008 s => 1250 Hz


 */
