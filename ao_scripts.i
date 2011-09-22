
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

func script_atmos_full(void)
{
  extern y_geom,y_tel,y_loop,y_atmos;
  extern g_atmos,g_target;
  extern ipupil;

  tabres = array(0.0f,4,4);

  filenames = [["atmos_1layer.par","atmos_1layer_vlt.par","atmos_1layer_20m.par","atmos_1layer_elt.par"],["atmos_1layer_5targets.par","atmos_1layer_vlt_5targets.par","atmos_1layer_20m_5targets.par","atmos_1layer_elt_5targets.par"],["atmos_4layers.par","atmos_4layers_vlt.par","atmos_4layers_20m.par","atmos_4layers_elt.par"],["atmos_12layers.par","atmos_12layers_vlt.par","atmos_12layers_20m.par","atmos_12layers_elt.par"]];
                 
  for (i=1;i<=4;i++) {
    for (j=1;j<=4;j++) {
      read_parfile,YOGA_AO_PARPATH+filenames(j,i);

      if (y_loop.niter == []) y_loop.niter = 1000;
      
      //here we just simulate atmosphere so pupdiam is fixed
      // arbitrarily.
      pupdiam = long(2*y_tel.diam*20.); // we assume subaps of 50cm and 20 phase pix per subaps

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
        /*
        if (cc % 10 == 0) {
          time_move = tac(mytime)/cc;
          //write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
        }
        */
      }
      write,"... Done !";
      // first dim is the telescope size
      // second dim is the type of parfile
      tabres(j,i) = tac(mytime)/y_loop.niter;
      
      //write,format="Average atmos move gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
    }
  }
  return tabres;
}



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

func script_wfs(filename,typeslp,verbose=)
{

  activeDevice,1;
  
  extern y_geom,y_tel,y_loop,y_atmos,y_wfs;
  extern g_atmos,g_target,g_wfs;
  extern ipupil;

  if (verbose == []) verbose = 1;
  
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

  if (verbose) write,"... Done with inits !";

  if (verbose) write,format="Starting loop on : %d iterations\n",y_loop.niter;
  
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
    if (verbose) { 
      if (cc % 10 == 0) {
        time_move = tac(mytime)/cc;
        write,format="\r Estimated remaining time : %.2f s",(y_loop.niter - cc)*time_move;
      }
    }
  }

  if (verbose) {
    write,"... Done !";
    write,format="Average wfs gpu time : %.4f s\n",tac(mytime)/y_loop.niter;
  }

  return tac(mytime)/y_loop.niter;
}

/*
// results
naos-like
diam=[4.,8.,20.,40.]
["1wfs8x8_1layer","1wfs1x16_1layer","1wfs40x40_1layer","1wfs80x80_1layer",]
4m          - 8m       - 20m      - 40m
naos1 = [0.000454419, 0.00117409, 0.0052497,  0.0253246]
naos2 = [0.000467373, 0.0011995,  0.00536267, 0.0257187]
naos3 = [0.000485916, 0.00128136, 0.00581861, 0.027374]
naos4 = [0.000470978, 0.00120419, 0.00538718, 0.0264424]
naos5 = [0.000496864, 0.00126634, 0.00575755, 0.0314704]
naos6 = [0.000473632, 0.001206,   0.00539444, 0.0263731]

naos-lgs-like
["1wfs8x8_1layer_lgs","1wfs1x16_1layer_lgs","1wfs40x40_1layer_lgs","1wfs80x80_1layer_lgs",]
4m          - 8m       - 20m     - 40m
naosLgs1 = [0.000763339, 0.00229221, 0.0121887, 0.218376]
naosLgs2 = [0.00077659 , 0.00231758, 0.0123006, 0.21854]
naosLgs3 = [0.00079519 , 0.0024001 , 0.0127567, 0.220232]
naosLgs4 = [0.000780166, 0.00232536, 0.01233  , 0.219253]
naosLgs5 = [0.000806483, 0.00238409, 0.0127017, 0.224476]
naosLgs6 = [0.000781823, 0.00232588, 0.0123366, 0.219242]

sphere-like
["sphere4m_1layer","sphere_1layer","sphere20m_1layer","sphere40m_1layer",]
4m          - 8m        - 20m      - 40m
sphere1 = [0.000343465, 0.000730022, 0.00327722, 0.0116978]
sphere2 = [0.000368606, 0.000797327, 0.0036923 , 0.0133116]
sphere3 = [0.000372487, 0.000863224, 0.00408113, 0.0148787]
sphere4 = [0.000374008, 0.000853045, 0.00405633, 0.014764]
sphere5 = [0.000432512, 0.0012298  , 0.00637453, 0.0240691]
sphere6 = [0.000377101, 0.000860718, 0.00405923, 0.0147827]


canary-like
["canary_1layer","canary8m_1layer","canary20m_1layer","canary40m_1layer",]
4m         - 8m       - 20m     - 40m
canary1 = [0.0021683 , 0.00573891, 0.0286964, 0.109454]
canary2 = [0.00222212, 0.00583664, 0.0291034, 0.110974]
canary3 = [0.00230654, 0.00608467, 0.0308139, 0.117758]
canary4 = [0.0022471 , 0.0059265 , 0.0298238, 0.113827]
canary5 = [0.00246509, 0.00669579, 0.03499  , 0.134702]
canary6 = [0.0022514 , 0.0059298 , 0.0297773, 0.113609]


*/

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
