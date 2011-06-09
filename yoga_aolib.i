plug_in,"yoga_ao";

require,"yoga.i";

/*
    _    ___              _                 _                 
   / \  / _ \    _____  _| |_ ___ _ __  ___(_) ___  _ __    _ 
  / _ \| | | |  / _ \ \/ / __/ _ \ '_ \/ __| |/ _ \| '_ \  (_)
 / ___ \ |_| | |  __/>  <| ||  __/ | | \__ \ | (_) | | | |  _ 
/_/   \_\___/   \___/_/\_\\__\___|_| |_|___/_|\___/|_| |_| (_)
                                                              
                                        
 _   _  ___   __ _  __ _     __ _  ___  
| | | |/ _ \ / _` |/ _` |   / _` |/ _ \ 
| |_| | (_) | (_| | (_| |  | (_| | (_) |
 \__, |\___/ \__, |\__,_|___\__,_|\___/ 
 |___/       |___/     |_____|          

*/

// atmosphere model
extern yoga_atmos;
extern init_tscreen;
extern get_tscreen;
extern get_spupil;
extern get_tscreen_config;
extern get_tscreen_update;
extern extrude_tscreen;

// targets
extern yoga_target;
extern target_addlayer;
extern target_atmostrace;
extern target_getphase;
extern target_getamplipup;

// phase
extern yoga_phase;
extern phase_copy;
extern phase_set;

// wfs
extern yoga_wfs;
extern wfs_initgs;
extern yoga_sensors;
extern sensors_initgs;
extern sensors_addlayer;
extern sensors_initarr;
extern sensors_compimg;
extern sensors_getimg;
extern sensors_getdata;

// global
extern move_atmos;
extern move_sky;
extern target_getimage;
