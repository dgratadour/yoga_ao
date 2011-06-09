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

//////////////////////////////////////////////////////////
// A basic set of display functions
//////////////////////////////////////////////////////////

func sky_win_init(pid1,pid2,pid3)
{
  extern sky_disp;

  sky_disp._xid=&[pid1,pid2,pid3];

  if (catch(0x08)) {
    sky_disp._gui_realized = 1;
  }
  
  if (!sky_disp._gui_realized) {
    stylename = yoga_ao_top+"/gs/yoga_ao.gs";
    window,(*sky_disp._wins)(1),dpi=sky_disp._defaultdpi,width=0,height=0,
      xpos=-2,ypos=-2,style=stylename,parent=pid1;
    limits,square=1;
    //palette,"gray.gp"; // need this if loadct is used!?
    /*
    widget_set_lut,sky_disp._lut,(*sky_disp._wins)(1),sky_disp._lut,
                   sky_disp._ncolors,sky_disp._itt,sky_disp._log_itt_dex,
                   sky_disp._rlut,sky_disp._glut,sky_disp._blut;
    */

    window,(*sky_disp._wins)(2),dpi=35,width=0,height=0,xpos=-2,ypos=-25,
      style=stylename,parent=pid2;
    limits,square=1;
    palette,"gray.gp"; 
   
    window,(*sky_disp._wins)(3),dpi=37,width=0,height=0,
      xpos=-2,ypos=-25,style=stylename,parent=pid3;
    limits,square=1;
    palette,"gray.gp"; 
   
  }

  sky_disp._gui_realized = 1;
}


/*
 _____                    _   
|_   _|_ _ _ __ __ _  ___| |_ 
  | |/ _` | '__/ _` |/ _ \ __|
  | | (_| | | | (_| |  __/ |_ 
  |_|\__,_|_|  \__, |\___|\__|
               |___/          
*/
func target_plot_update(ntarget)
{
  if (y_target == []) return;

  window,(*sky_disp._wins)(3);fma;
  plmk,*y_target.ypos,*y_target.xpos,marker=5;
  plmk,(*y_target.ypos)(ntarget),(*y_target.xpos)(ntarget),marker=5,color="red";
  limits,-40,40;
  range,-40,40;
  xytitles,"xpos (arcsec)","ypos (arcsec)",[0.02,0.02];
  yoga_pltitle,"Asterism",[0.,-0.015];
  
  
}
func create_target(ntargets)
{
  extern y_target;
  
  y_target  = target_struct();
  
  y_target.ntargets = ntargets;
  if (ntargets > 1) {
    y_target.xpos     = &(random_n(ntargets));
    y_target.ypos     = &(random_n(ntargets));
  } else {
    y_target.xpos     = &([0.]);
    y_target.ypos     = &([0.]);
  }
  y_target.mag      = &(array(5.,ntargets));
  y_target.lambda   = &(array(1.6,ntargets));

  for (cc=1;cc<=100;cc++)
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').remove_text(%d)",0);
  for (cc=1;cc<=ntargets;cc++)
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",0,swrite(format="Source # %d",ntargets-cc+1));
  
  update_target_prop,0;
  pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
}

func update_target_prop(ntarget)
{
  extern sky_disp,y_target;
  
  if (y_target == []) {
    pyk,swrite(format=sky_disp._cmd+"y_init_target(%d)",1);
  } else {
    if (ntarget > y_target.ntargets) error,"not a valid target";
    else {
      pyk,swrite(format=sky_disp._cmd+"y_update_target_gui(%f, %f, %f, %f)",
                 (*y_target.xpos)(ntarget),(*y_target.ypos)(ntarget),
                 (*y_target.lambda)(ntarget),(*y_target.mag)(ntarget));
      pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
    }
  }
}

func update_ntargets(ntargets)
{
  extern y_target;
  
  if (y_target != []) {
    if (ntargets > y_target.ntargets) {
      y_target.ntargets += 1;
      y_target.xpos      = &(_((*y_target.xpos),0.));
      y_target.ypos      = &(_((*y_target.ypos),0.));
      y_target.lambda    = &(_((*y_target.lambda),1.6));
      y_target.mag       = &(_((*y_target.mag),5));
      
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",
                 ntargets-1,swrite(format="Source # %d",y_target.ntargets));
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",y_target.ntargets-1);
    
    } else if (ntargets <  y_target.ntargets) {
      y_target.ntargets -= 1;
      y_target.xpos      = &((*y_target.xpos)(1:-1));
      y_target.ypos      = &((*y_target.ypos)(1:-1));
      y_target.lambda    = &((*y_target.lambda)(1:-1));
      y_target.mag       = &((*y_target.mag)(1:-1));

      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').remove_text(%d)",y_target.ntargets);
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",y_target.ntargets-1);
    }
    pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
    if (ntargets != y_target.ntargets) update_ntargets,ntargets;
  }
}

func init_target_prop(ntarget,mag,xpos,ypos,lambda)
{
  extern y_target;

  if (y_target == []) return;
  else {
    (*y_target.mag)(ntarget) = mag;
    (*y_target.xpos)(ntarget) = xpos;
    (*y_target.ypos)(ntarget) = ypos;
    (*y_target.lambda)(ntarget) = lambda;
    pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
  }
}

func remove_target(ntarget)
{
  extern y_target;

  if (y_target == []) return;
  else {
    if (y_target.ntargets > 1) {
      //      y_atmos.nscreens -=1;      
      if (ntarget != y_target.ntargets) {
        if (ntarget == 1) {
          *y_target.xpos = roll(*y_target.xpos,-1);
          *y_target.ypos = roll(*y_target.ypos,-1);
          *y_target.lambda = roll(*y_target.lambda,-1);
          *y_target.mag = roll(*y_target.mag,-1);
        } else {
          xpos = *y_target.xpos;
          y_target.xpos = &(_(xpos(1:ntarget-1),xpos(ntarget+1:),xpos(ntarget)));
          ypos = *y_target.ypos;
          y_target.ypos = &(_(ypos(1:ntarget-1),ypos(ntarget+1:),ypos(ntarget)));
          lambda = *y_target.lambda;
          y_target.lambda = &(_(lambda(1:ntarget-1),lambda(ntarget+1:),lambda(ntarget)));
          mag = *y_target.mag;
          y_target.mag = &(_(mag(1:ntarget-1),mag(ntarget+1:),mag(ntarget)));
        }
      }
      //update_nlayers(pupdiam,y_atmos.nscreens-1); 
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('ntargets').set_value(%d)",y_target.ntargets-1);
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('layer_select').set_active(%d)",ntarget-1);
      //pyk,swrite(format=sky_disp._cmd+"y_layers_plot_update(%d)",1);
    }
  }
}

func load_default_target(tconf)
{
  extern y_target;
  
  y_target  = target_struct(); // clean start

  if (tconf == 1) { // one source ...  
    y_target.ntargets = 1;
    y_target.xpos     = &([0.0]);
    y_target.ypos     = &([0.0]);
    y_target.lambda   = &([1.6]);
    y_target.mag      = &([5.]);
  }
  if (tconf == 2) { // two sources ...  
    y_target.ntargets = 2;
    y_target.xpos     = &([20.0,-20.0]);
    y_target.ypos     = &([0.0,0.0]);
    y_target.lambda   = &([1.6,1.6]);
    y_target.mag      = &([5.,5.]);
  }
  if (tconf == 3) { // three sources ...  
    y_target.ntargets = 3;
    y_target.xpos     = &([20.0,-10.0,-10.0]);
    y_target.ypos     = &([0.0,17.3,-17.3]);
    y_target.lambda   = &([1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.]);
  }
  if (tconf == 4) { // four sources ...  
    y_target.ntargets = 4;
    y_target.xpos     = &([20.0,0.0,-20.0,0.0]);
    y_target.ypos     = &([0.0,20.0,0.0,-20.0]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.]);
  }
  if (tconf == 5) { // five sources ...  
    y_target.ntargets = 5;
    y_target.xpos     = &([20.0,6.2,-16.2,-16.2,6.2]);
    y_target.ypos     = &([0.0,19.0,11.7,-11.7,-19.0]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.,5.]);
  }
  if (tconf == 6) { // six sources ...  
    y_target.ntargets = 6;
    y_target.xpos     = &([20.0,10.0,-10.0,-20.0,-10.,10.]);
    y_target.ypos     = &([0.0,17.3,17.3,0.0,-17.3,-17.3]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.,5.,5.]);
  }
  if (tconf == 7) { // seven sources ...  
    y_target.ntargets = 7;
    y_target.xpos     = &([20.0,10.,-10.,-20.0,-10.,10.,0.0]);
    y_target.ypos     = &([0.0,17.3,17.3,0.0,-17.3,-17.3,0.0]);
    y_target.lambda   = &([1.6,1.6,1.6,1.6,1.6,1.6,1.6]);
    y_target.mag      = &([5.,5.,5.,5.,5.,5.,5.]);
  }
  if (tconf < 8) {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').remove_text(%d)",0);
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').insert_text(%d,'%s')",0,swrite(format="Source # %d",y_target.ntargets-cc+1));
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('target_select').set_active(%d)",0);
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('ntargets').set_value(%d)",y_target.ntargets);
    
    pyk,swrite(format=sky_disp._cmd+"y_target_plot_update(%d)",1);
  }
}

func init_all(pupdiam,zenith,teldiam,cobs,r0,freq)
{
  extern y_atmos,y_target;
  extern g_target;
  
  init_gatmos,pupdiam,zenith,teldiam,cobs,r0/100.,freq;

  target_init;
  
  update_main_display1,["atmos","target","image"];
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
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').remove_text(%d)",0);
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",0,
               "Phase - Atmos");
  }
  if (anyof(type == "target")) {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",1,
               "Phase - Target");
  }
  if (anyof(type == "image")) {
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').insert_text(%d,'%s')",2,
               "Image - Target");
  }
  
  pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_type').set_active(%d)",0);
}

func update_main_display2(type)
{
  extern y_atmos;

  if (type == "Phase - Atmos") {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').remove_text(%d)",0);
    for (cc=1;cc<=y_atmos.nscreens;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Layer # %d",y_atmos.nscreens-cc+1));  
  }
  if (type == "Phase - Target") {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').remove_text(%d)",0);
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));  
  }
  if (type == "Image - Target") {
    for (cc=1;cc<=100;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').remove_text(%d)",0);
    for (cc=1;cc<=y_target.ntargets;cc++)
      pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').insert_text(%d,'%s')",0,
                 swrite(format="Source # %d",y_target.ntargets-cc+1));  
  }
  pyk,swrite(format=sky_disp._cmd+"glade.get_widget('winselect_number').set_active(%d)",0);
}

func update_main(type,nlayer)
{
  extern y_atmos,g_atmos,sky_disp;
  if (nlayer < 0) return;
  
  if (type == "Phase - Atmos") {
    mscreen = get_tscreen(g_atmos,(*y_atmos.alt)(nlayer+1));
    window,(*sky_disp._wins)(1);fma;limits;
    pli,mscreen;
  }
  if (type == "Phase - Target") {
    target_atmostrace,g_target,nlayer,g_atmos;
    mscreen = target_getphase(g_target,nlayer);
    window,(*sky_disp._wins)(1);fma;limits;
    pli,mscreen;
  }
  if (type == "Image - Target") {
    mscreen = target_getimage(g_target,g_atmos,nlayer);
    window,(*sky_disp._wins)(1);fma;limits;
    pli,eclat(mscreen);
  }
}

func start_sky_loop
{
  extern skyiter,skyloop;

  skyiter = 0;
  skyloop = 1;
  animate,1;
  sky_loop;
}

func sky_loop(one)
{
  extern skyloop,skydisp_type,skydisp_num;
  extern g_atmos,g_target;
  extern y_atmos;
  extern time_move;

  if (skyloop) {
    if (g_atmos == []) return;

    mytime = tic();

    //if (g_target == []) extrude_tscreen,g_atmos,(*y_atmos.alt)(1);
    if (g_target == []) move_atmos,g_atmos;
    else move_sky,g_atmos,g_target;

    
    if (skydisp_type != [])
      update_main,skydisp_type,skydisp_num;    
    
    if (skyiter < 1) time_move = tac(mytime);
    else {
      time_move += 0.01*tac(mytime);
      time_move /= 1.01;
    }

    skyiter ++;

    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('progressbar_sky').set_fraction(%f)",float((skyiter%1000)/1000.));
    mtext = swrite(format="Framerate : %.2f",1./time_move);
    pyk,swrite(format=sky_disp._cmd+"glade.get_widget('progressbar_sky').set_text('%s')",mtext);
    
    if (!one) after,0.001,sky_loop;
  } else return;
}


//////////////////////////////////////////////////////////
//              **********************                 //
////////         END OF ROUTINES DEFINITION      ////////
//              **********************                 //
//////////////////////////////////////////////////////////

// start standalone version if called from shell  
//pyk_debug=1;

arg_sky = get_argv();
sky_disp = display_struct();

sky_disp._cmd = "sky.";

if (anyof(strmatch(arg_sky,"widget_sky.i")) || get_env("EMACS")=="t" ) {
  sky_disp._cmd = "";
  python_exec = yoga_ao_top+"/widgets/widget_sky.py";
  pyk_cmd=[python_exec];
  if (!_pyk_proc) _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  write,"widget_sky  ready";
  write,"standalone version";
}

sky_disp._wins            = &[10,11,12];
sky_disp._defaultdpi      = 85;    // change size of spydr graphic area
sky_disp._ncolors         = 200;
sky_disp._lut             = 0;     // default LUT index [0-41]
sky_disp._xytitles_adjust = &[0.012,0.019]; // X and Y notch axis titles in main area
sky_disp._invertlut       = 0;
pldefault,opaque=1,marks  =1;
sky_disp._gui_realized    = 0;
sky_disp._zoom            = 1;

skydisp_type = "";
skydisp_num  = 0;
skyloop      = 0;
skyiter      = 0;

atmos_disp = sky_disp;
