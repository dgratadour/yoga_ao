/*
*/
// Environment check
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;
/*
plug_path = yoga_ao_top;
grow,plug_path,"/home/brujo/yorick/yoga/trunk";
plug_dir,plug_path;
*/
// load necessary files : yorick-python wrapper and styc_utils
require,yoga_ao_top+"/ywidgets/pyk.i";
require,yoga_ao_top+"/ywidgets/widgets_utils.i";
require,yoga_ao_top+"/ywidgets/atmos_utils.i";
require,yoga_ao_top+"/yoga_ao.i";


func wfs_win_init(pid1,pid2)
{
  extern wfs_disp;

  wfs_disp._xid=&[pid1,pid2];

  if (catch(0x08)) {
    wfs_disp._gui_realized = 1;
  }
  
  if (!wfs_disp._gui_realized) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs";
    window,(*wfs_disp._wins)(1),dpi=wfs_disp._defaultdpi,width=0,height=0,
      xpos=-4,ypos=-4,style=stylename,parent=pid1;
    limits,square=1;

    window,(*wfs_disp._wins)(2),dpi=35,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid2;
    limits,square=1;
    palette,"gray.gp"; 
   
  }

  wfs_disp._gui_realized = 1;
}

func update_wfs_prop(numwfs)
{
  extern wfs_disp,y_wfs;
  
  if (y_wfs == []) {
    pyk,swrite(format=sky_disp._cmd+"y_init_wfs(%d)",1);
  } else {
    if (numwfs > numberof(y_wfs)) error,"not a valid target";
    else {
      pyk,swrite(format=wfs_disp._cmd+"y_update_wfs_gui(%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f)",
                 (y_wfs.nxsub)(numwfs),(y_wfs.npix)(numwfs),
                 (y_wfs.pixsize)(numwfs),(y_wfs.fracsub)(numwfs),
                 (y_wfs.xpos)(numwfs),(y_wfs.ypos)(numwfs),(y_wfs.lambda)(numwfs),
                 (y_wfs.gsmag)(numwfs),log10((y_wfs.zerop)(numwfs)),
                 (y_wfs.optthroughput)(numwfs),(y_wfs.noise)(numwfs));
    }
  }
}

func init_wfs_prop(numwfs,nsub,npix,pixsize,mag,xpos,ypos,lambda,frac,zp,throughput,noise)
{
  extern y_wfs;

  if (y_wfs == []) return;
  else {
    y_wfs(numwfs).nxsub         = nsub;
    y_wfs(numwfs).npix          = npix;
    y_wfs(numwfs).pixsize       = pixsize;
    y_wfs(numwfs).fracsub       = frac;
    y_wfs(numwfs).xpos          = xpos;
    y_wfs(numwfs).ypos          = ypos;
    y_wfs(numwfs).lambda        = lambda;
    y_wfs(numwfs).gsmag         = mag;
    y_wfs(numwfs).optthroughput = throughput;
    y_wfs(numwfsxs).zerop       = 10^zp;
    y_wfs(numwfsxs).noise       = noise;
  }
}

func init_wfs_prop_lgs(gsalt,lltx,llty,power,wreturn,proftype,beam)
{
  extern y_wfs;

  if (y_wfs == []) return;
  else {
    y_wfs(numwfs).gsalt            = gsalt*1.e3;
    y_wfs(numwfs).lltx             = lltx;
    y_wfs(numwfs).llty             = llty;
    y_wfs(numwfs).laserpower       = power;
    y_wfs(numwfs).lgsreturnperwatt = wreturn;
    y_wfs(numwfs).proftype         = proftype;
    y_wfs(numwfs).beamsize         = beam;
  }
}


func create_wfs(numwfs)
{
  extern y_wfs;
  
  y_wfs  = array(wfs_struct(),numwfs);
  
  if (ntargets > 1) {
    xpos     = random_n(numwfs);
    ypos     = random_n(numwfs);
  } else {
    xpos     = 0.;
    ypos     = 0.;
  }
  for (i=1;i<=numwfs;i++) {
      y_wfs(i).nxsub         = 7;
      y_wfs(i).npix          = 7;
      y_wfs(i).pixsize       = 0.5;
      y_wfs(i).fracsub       = 0.8;
      y_wfs(i).xpos          = xpos(i);
      y_wfs(i).ypos          = ypos(i);
      y_wfs(i).lambda        = 0.5;
      y_wfs(i).gsmag         = 5.;
      y_wfs(i).optthroughput = 0.5;
      y_wfs(i).zerop         = 1.e11;
      y_wfs(i).noise         = -1;
  }
  
  for (cc=1;cc<=100;cc++)
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').remove_text(%d)",0);
  for (cc=1;cc<=ntargets;cc++)
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",0,swrite(format="WFS # %d",numwfs-cc+1));
  
  update_wfs_prop,0;
}

func update_nwfs(numwfs)
{
  extern y_wfs;
  
  if (y_wfs != []) {
    if (numwfs > numberof(y_wfs)) {
      tmp = y_wfs;
      tmp2 = array(wfs_struct,numberof(y_wfs)+1);
      tmp2(1:numberof(y_wfs)) = y_wfs;
      y_wfs = tmp2;
      
      y_wfs(0).nxsub         = 7;
      y_wfs(0).npix          = 7;
      y_wfs(0).pixsize       = 0.5;
      y_wfs(0).fracsub       = 0.8;
      y_wfs(0).xpos          = 0.0;
      y_wfs(0).ypos          = 0.0;
      y_wfs(0).lambda        = 0.5;
      y_wfs(0).gsmag           = 5.;
      y_wfs(0).optthroughput = 0.5;
      y_wfs(0).zerop         = 1.e11;
      y_wfs(0).noise         = -1;
      
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",
                 numberof(y_wfs)-1,swrite(format="WFS # %d",numberof(y_wfs)));
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",numberof(y_wfs)-1);
    
    } else if (numwfs <  numberof(y_wfs)) {
      y_wfs = y_wfs(1:-1);

      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').remove_text(%d)",numberof(y_wfs));
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",numberof(y_wfs)-1);
    }
     if (numwfs != numberof(y_wfs)) update_nwfs,numwfs;
  }
}

func remove_wfs(numwfs)
{
  extern y_wfs;

  if (y_wfs == []) return;
  else {
    if (numberof(y_wfs)>1) {
      if (numwfs != numberof(y_wfs)) {
        if (numwfs == 1) {
          tmp = y_wfs;
          tmp(1:numberof(y_wfs)-1) = y_wfs(2:);
          y_wfs = tmp;
        } else {
          tmp = y_wfs;
          tmp(1:numwfs-1) = y_wfs(1:numwfs-1);
          tmp(numwfs:-1) = y_wfs(numwfs+1:);
          y_wfs = tmp;
        }
      }
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs)-1);
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",numberof(y_wfs)-1);
    }
  }
}


func load_default_wfs(tconf)
{
  extern y_wfs;
  
  if (tconf == 1) { // one wfs ...  
    y_wfs          = array(wfs_struct(),1); // clean start
    y_wfs(1).nxsub         = 7;
    y_wfs(1).npix          = 7;
    y_wfs(1).pixsize       = 0.5;
    y_wfs(1).fracsub       = 0.8;
    y_wfs(1).xpos          = 0.0;
    y_wfs(1).ypos          = 0.0;
    y_wfs(1).lambda        = 0.5;
    y_wfs(1).gsmag         = 5.;
    y_wfs(1).optthroughput = 0.5;
    y_wfs(1).zerop         = 1.e11;
    y_wfs(1).noise         = -1;
  }
  if (tconf < 8) {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').remove_text(%d)",0);
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').insert_text(%d,'%s')",0,swrite(format="WFS # %d",numberof(y_wfs)-cc+1));
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('wfs_select').set_active(%d)",0);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('nwfs').set_value(%d)",numberof(y_wfs));
    
  }
}


func start_wfs_loop
{
  extern wfsiter,wfsloop;

  wfsiter = 0;
  wfsloop = 1;
  animate,1;
  wfs_loop;
}

func wfs_loop(one)
{
  extern wfsloop,wfsdisp_type,wfsdisp_num;
  extern g_atmos;
  extern y_atmos;
  extern time_move;

  if (wfsloop) {
    if (g_atmos == []) return;

    mytime = tic();

    //if (g_target == []) extrude_tscreen,g_atmos,(*y_atmos.alt)(1);
    move_atmos,g_atmos;
    if ((y_wfs != []) && (g_wfs != [])) {
      for (i=1;i<=numberof(y_wfs);i++) {
        sensors_trace,g_wfs,i-1,g_atmos;
        sensors_compimg,g_wfs,i-1;
      }
    }
    
    if (wfsdisp_type != [])
      update_main,wfsdisp_type,wfsdisp_num;    
    
    if (wfsiter < 1) time_move = tac(mytime);
    else {
      time_move += 0.01*tac(mytime);
      time_move /= 1.01;
    }

    wfsiter ++;

    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('progressbar_wfs').set_fraction(%f)",float((wfsiter%1000)/1000.));
    mtext = swrite(format="Framerate : %.2f",1./time_move);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('progressbar_wfs').set_text('%s')",mtext);
    
    if (!one) after,0.001,wfs_loop;
  } else return;
}

func init_all(zenith,teldiam,cobs,r0,freq)
{
  extern y_atmos,y_wfs;
  extern g_wfs;

  y_atmos.r0 = r0;
  if (y_tel == []) y_tel   = tel_struct();
  if (y_loop == []) y_loop = loop_struct();
  
  y_tel.diam         = teldiam;
  y_tel.cobs         = cobs;
  y_loop.ittime      = 1.0/freq;
  
  wfs_init;
  
  init_gatmos,y_geom.pupdiam,zenith,teldiam,cobs,r0,freq;

  update_main_display1,["atmos","wfs","image"];
  
  target_init;

  y_target = [];
}
/*
 __  __       _                                _ 
|  \/  | __ _(_)_ __    _ __   __ _ _ __   ___| |
| |\/| |/ _` | | '_ \  | '_ \ / _` | '_ \ / _ \ |
| |  | | (_| | | | | | | |_) | (_| | | | |  __/ |
|_|  |_|\__,_|_|_| |_| | .__/ \__,_|_| |_|\___|_|
                       |_|                       
*/
func update_main_display1(type)
{
  if (anyof(type == "atmos")) {
    for (cc=1;cc<=4;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').remove_text(%d)",0);
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",0,
               "Phase - Atmos");
  }
  
  if (anyof(type == "wfs")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",1,
               "Phase - WFS");
  }
  
  if (anyof(type == "image")) {
    pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",2,
               "Image - WFS");
  }
  
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_type').set_active(%d)",0);
}

func update_main_display2(type)
{
  extern y_atmos;

  if (type == "Phase - Atmos") {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').remove_text(%d)",0);
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",y_atmos.nscreens-cc+1));  
  }
  
  if (type == "Phase - WFS") {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').remove_text(%d)",0);
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  if (type == "Image - WFS") {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').remove_text(%d)",0);
    for (cc=1;cc<=numberof(y_wfs);cc++)
      pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="WFS # %d",numberof(y_wfs)-cc+1));  
  }
  
  pyk,swrite(format=wfs_disp._cmd+"glade.get_widget('winselect_number').set_active(%d)",0);
}

func update_main(type,nlayer)
{
  extern y_atmos,g_atmos,wfs_disp;
  if (nlayer < 0) return;
  
  if (type == "Phase - Atmos") {
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(nlayer+1));
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
  }
  
  if (type == "Phase - WFS") {
    sensors_trace,g_wfs,nlayer,g_atmos;
    //sensors_compimg,g_wfs,nlayer;
    mscreen = sensors_getdata(g_wfs,nlayer,"phase");
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
  }
  
  if (type == "Image - WFS") {
    sensors_trace,g_wfs,nlayer,g_atmos;
    sensors_compimg,g_wfs,nlayer;
    mscreen = sensors_getimg(g_wfs,nlayer);
    window,(*wfs_disp._wins)(1);fma;
    pli,mscreen;
  }
}

//////////////////////////////////////////////////////////
//              **********************                 //
////////         END OF ROUTINES DEFINITION      ////////
//              **********************                 //
//////////////////////////////////////////////////////////

// start standalone version if called from shell  
//pyk_debug=1;

arg_wfs = get_argv();
wfs_disp = display_struct();

wfs_disp._cmd = "wfs.";

if (anyof(strmatch(arg_wfs,"widget_wfs.i")) || get_env("EMACS")=="t" ) {
  wfs_disp._cmd = "";
  python_exec = yoga_ao_top+"/widgets/widget_wfs.py";
  pyk_cmd=[python_exec];
  if (!_pyk_proc) _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  write,"widget_wfs  ready";
  write,"standalone version";
}

wfs_disp._wins            = &[13,14];
wfs_disp._defaultdpi      = 145;    // change size of spydr graphic area
wfs_disp._ncolors         = 200;
wfs_disp._lut             = 0;     // default LUT index [0-41]
wfs_disp._xytitles_adjust = &[0.012,0.019]; // X and Y notch axis titles in main area
wfs_disp._invertlut       = 0;
pldefault,opaque=1,marks  =1;
wfs_disp._gui_realized    = 0;
wfs_disp._zoom            = 1;

wfsdisp_type = "";
wfsdisp_num  = 0;
wfsloop      = 0;
wfsiter      = 0;

atmos_disp = wfs_disp;
