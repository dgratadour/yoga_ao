#include <yoga_turbu.h>
#include <yoga_target.h>
#include <yoga_phase.h>
#include <yoga_wfs.h>
#include <yapi.h>

int move_atmos(yoga_atmos *atmos,yoga_target *target) {

  map<float, yoga_tscreen *>::iterator p;
  map<type_screen,yoga_obj<float> *>::iterator pp;
  p = atmos->d_screens.begin();
  yoga_tscreen *tmp;
  int cpt=1;
  while (p != atmos->d_screens.end()) {
    tmp = p->second;
    tmp->accumx += tmp->deltax;
    tmp->accumy += tmp->deltay;

    int deltax = (int)tmp->accumx;
    int deltay = (int)tmp->accumy;

    for (int cc=0;cc<deltax;cc++) tmp->extrude(1);
    tmp->accumx -= deltax;
    for (int dd=0;dd<target->ntargets;dd++) {
      pp = target->d_targets.at(dd)->rays_posx.find(make_pair("atmos",tmp->altitude));
      if (pp != target->d_targets.at(dd)->rays_posx.end()) {
	yoga_plus(pp->second->d_data,(tmp->deltax-deltax),pp->second->nb_elem);
      }
    }

    for (int cc=0;cc<deltay;cc++) tmp->extrude(0);
    tmp->accumy -= deltay;
    for (int dd=0;dd<target->ntargets;dd++) {
      pp = target->d_targets.at(dd)->rays_posy.find(make_pair("atmos",tmp->altitude));
      if (pp != target->d_targets.at(dd)->rays_posy.end()) {
	yoga_plus(pp->second->d_data,(tmp->deltay-deltay),pp->second->nb_elem);
      }
    }
    p++;cpt++;
  }

}

int move_atmos(yoga_atmos *atmos) {

  map<float, yoga_tscreen *>::iterator p;
  p = atmos->d_screens.begin();
  yoga_tscreen *tmp;
  while (p != atmos->d_screens.end()) {
    p->second->accumx += p->second->deltax;
    p->second->accumy += p->second->deltay;

    int deltax = (int)p->second->accumx;
    int deltay = (int)p->second->accumy;
    for (int cc=0;cc<deltax;cc++) p->second->extrude(1);
    p->second->accumx -= deltax;
    for (int cc=0;cc<deltay;cc++) p->second->extrude(0);
    p->second->accumy -= deltay;
    p++;
  }

}

/*
__   __         _      _         _    ____ ___ 
\ \ / /__  _ __(_) ___| | __    / \  |  _ \_ _|
 \ V / _ \| '__| |/ __| |/ /   / _ \ | |_) | | 
  | | (_) | |  | | (__|   <   / ___ \|  __/| | 
  |_|\___/|_|  |_|\___|_|\_\ /_/   \_\_|  |___|
                                               
 */


extern "C" {

  /*
 _   ____                           
| |_/ ___|  ___ _ __ ___  ___ _ __  
| __\___ \ / __| '__/ _ \/ _ \ '_ \ 
| |_ ___) | (__| | |  __/  __/ | | |
 \__|____/ \___|_|  \___|\___|_| |_|
                                    
   */

  typedef struct tscreen_struct {
    void *yoga_tscreen;
    int device;    
  } tscreen_struct;


  void tscreen_free(void *obj) {
    tscreen_struct *handler = (tscreen_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_tscreen *tscreen_obj_handler = (yoga_tscreen *)(handler->yoga_tscreen);
      delete tscreen_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void tscreen_print(void *obj)
  {
    tscreen_struct *handler = (tscreen_struct *)obj;
    yoga_tscreen *tscreen_obj_handler = (yoga_tscreen *)(handler->yoga_tscreen);
    cout << "Yoga tScreen Object : " << endl;
    cout << "Screen @ alt. " << tscreen_obj_handler->altitude
	 << " m |" << " wind speed : " << tscreen_obj_handler->windspeed << " m/s |"
	 << " wind dir. : " << tscreen_obj_handler->winddir << " deg |"
	 << "r_0 : "<<  pow(tscreen_obj_handler->amplitude,-6./5.) << " m"<<endl;
  }


  string y_tscreenName("yTscreen Object");
  static y_userobj_t yTscreen =
    {const_cast<char*>(y_tscreenName.c_str()), 
     &tscreen_free, 0, 0, 0, 0};

  /*
    _   _                       
   / \ | |_ _ __ ___   ___  ___ 
  / _ \| __| '_ ` _ \ / _ \/ __|
 / ___ \ |_| | | | | | (_) \__ \
/_/   \_\__|_| |_| |_|\___/|___/
                                
   */

  typedef struct atmos_struct {
    void *yoga_atmos;
    int device;    
  } atmos_struct;

  void atmos_free(void *obj) {
    atmos_struct *handler = (atmos_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_atmos *atmos_obj_handler = (yoga_atmos *)(handler->yoga_atmos);
      delete atmos_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void atmos_print(void *obj)
  {
    atmos_struct *handler = (atmos_struct *)obj;
    yoga_atmos *atmos_obj_handler = (yoga_atmos *)(handler->yoga_atmos);
    map<float, yoga_tscreen *>::iterator p;
    p = atmos_obj_handler->d_screens.begin();
    cout << "Yoga Atmos Object : " << endl;
    cout << "Contains " << atmos_obj_handler->nscreens << " turbulent screen(s) : " << endl;
    cout << "r_0 = " << atmos_obj_handler->r0 << " m" << endl;
    int i =0;
    yoga_tscreen *tmp;
    while (p != atmos_obj_handler->d_screens.end()) {
      tmp = p->second;
      cout << "Screen # " << i+1 << " @ alt. " << tmp->altitude
	   << " m |" << " wind speed : " << tmp->windspeed << " m/s |"
	   << " wind dir. : " << tmp->winddir << " deg |"
	   << "r_0 : "<<  pow(tmp->amplitude,-6/5.) << " m"
	   << " deltax : " << tmp->deltax << " pix | deltay : "
	   << tmp->deltay << " pix "<< endl;
      p++; i++;
    }
    //delete tmp;
  }


  string y_AtmosName("yAtmos Object");
  static y_userobj_t yAtmos =
    {const_cast<char*>(y_AtmosName.c_str()), 
     &atmos_free, &atmos_print, 0, 0, 0};

  void
  Y_yoga_atmos(int argc)
  {    
    long ntot;
    long dims[Y_DIMSIZE];
    try { 
      //(int nscreens,T_screen r0,long *size,float *frac, float *altitude, float *windspeed, float *winddir)
      int odevice = activeDevice;
      int nscreens = ygets_i(argc-1);
      float *r0 = ygeta_f(argc-2, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens size");

      long *size;
      size = ygeta_l(argc-3, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens size");

      long *size2;
      size2 = ygeta_l(argc-4, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens size");

      float *frac;
      frac = ygeta_f(argc-5, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens frac");

      float *alt;
      alt = ygeta_f(argc-6, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens alt");

      float *wspeed;
      wspeed = ygeta_f(argc-7, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens speed");

      float *wdir;
      wdir = ygeta_f(argc-8, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens speed");

      float *deltax;
      deltax = ygeta_f(argc-9, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens speed");

       float *deltay;
      deltay = ygeta_f(argc-10, &ntot,dims);
      if (ntot != nscreens) y_error("wrong dimension for screens speed");

      float *pupil;
      pupil = ygeta_f(argc-11, &ntot,dims);
      //if (ntot != nscreens) y_error("wrong dimension for screens speed");

     if (argc > 11) odevice = ygets_i(argc-12);
      if (odevice != activeDevice) {
	if (odevice > _nbDevice()-1) {
	  cout << "### Invalid Device Id : " << odevice <<" Your system has only " << 
	    _nbDevice() << " device(s) available" << endl;
	  cout << "Using active Device with Id : " << activeDevice << endl;
	} else {
	  cudaDeviceProp deviceProp;
	  cutilSafeCall( cudaGetDeviceProperties(&deviceProp, odevice) );
	  yoga_setDevice(odevice);
	}
      }

      atmos_struct *handle=(atmos_struct *)ypush_obj(&yAtmos, sizeof(atmos_struct));
      handle->device = odevice;
      
      handle->yoga_atmos = new yoga_atmos(nscreens,(float *)r0,(long *)size,(long *)size2,(float *)alt,(float *)wspeed,(float *)wdir,(float *)deltax,(float *)deltay,(float *)pupil);
      
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    } catch( ... ) { 
      stringstream buf;
      buf << "unknown error with yoga_atmos construction in "<<__FILE__ << "@" << __LINE__ << endl;
      y_error(buf.str().c_str()); 
    }
  }
 
  void
  Y_init_tscreen(int argc)
  // here we init the atmos structure
  // args are : the yoga_atmos, nscreen A, B, istencilx,istencily and seed
  {    
    long ntot;
    long dims[Y_DIMSIZE];
    atmos_struct *handler = (atmos_struct *)yget_obj(argc-1,&yAtmos);
    yoga_atmos *atmos_obj_handler = (yoga_atmos *)(handler->yoga_atmos);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    float altitude = ygets_f(argc-2);

    float *a = ygeta_f(argc-3, &ntot,dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_A->nb_elem)
      y_error("wrong size for A");
    float *b = ygeta_f(argc-4, &ntot,dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_B->nb_elem)
      y_error("wrong size for B");
    unsigned int *istencilx = (unsigned int *)ygeta_i(argc-5, &ntot,dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_istencilx->nb_elem)
      y_error("wrong size for istencilx");
    unsigned int *istencily = (unsigned int *)ygeta_i(argc-6, &ntot,dims);
    if (ntot != atmos_obj_handler->d_screens[altitude]->d_istencily->nb_elem)
      y_error("wrong size for istencily");

    int seed = ygets_i(argc-7);

    atmos_obj_handler->init_screen(altitude,a,b,istencilx,istencily,seed);
  }

  void
  Y_get_tscreen(int argc)
  { 

    if (yarg_subroutine()) y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *)yget_obj(argc-1,&yAtmos);
      if (handle->device != activeDevice) yoga_setDevice(handle->device);
      float alt = ygets_f(argc-2);
      yoga_atmos *atmos_handler = (yoga_atmos *)handle->yoga_atmos;
      yObjS *yoga_obj_handler = (yObjS *)(atmos_handler->d_screens[alt]->d_tscreen->d_screen);
      float *data = ypush_f(yoga_obj_handler->dims_data);
      yoga_obj_handler->device2host(data);
    }
  }

  void
  Y_get_spupil(int argc)
  { 

    if (yarg_subroutine()) y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *)yget_obj(argc-1,&yAtmos);
      if (handle->device != activeDevice) yoga_setDevice(handle->device);
      yoga_atmos *atmos_handler = (yoga_atmos *)handle->yoga_atmos;
      if (atmos_handler->d_screens.find(0.0f) != atmos_handler->d_screens.end()) {
	yObjS *yoga_obj_handler = (yObjS *)(atmos_handler->d_pupil);
	float *data = ypush_f(yoga_obj_handler->dims_data);
	yoga_obj_handler->device2host(data);
      }
    }
  }

  void
  Y_extrude_tscreen(int argc)
  { 

    if (yarg_subroutine())  {
      atmos_struct *handle = (atmos_struct *)yget_obj(argc-1,&yAtmos);
      if (handle->device != activeDevice) yoga_setDevice(handle->device);
      float alt = ygets_f(argc-2);
      yoga_atmos *atmos_handler = (yoga_atmos *)handle->yoga_atmos;
      yoga_tscreen *tscreen_handler = (yoga_tscreen *)(atmos_handler->d_screens.find(alt)->second);
      int dir;
      if (argc>2) dir = ygets_i(argc-3);
      else dir = 1;
     tscreen_handler->extrude(dir);
   } else y_error("can only be called as a subroutine");
  }

  void
  Y_get_tscreen_config(int argc)
  { 

    if (yarg_subroutine())y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *)yget_obj(argc-1,&yAtmos);
      if (handle->device != activeDevice) yoga_setDevice(handle->device);
      float alt = ygets_f(argc-2);
      yoga_atmos *atmos_handler = (yoga_atmos *)handle->yoga_atmos;
      yoga_tscreen *tscreen_handler = (yoga_tscreen *)(atmos_handler->d_screens[alt]);
      char *type_data = ygets_q(argc-3);
      if (strcmp(type_data, "A")==0) {
	float *data = ypush_f(tscreen_handler->d_A->dims_data);
	tscreen_handler->d_A->device2host(data);
      } else if (strcmp(type_data, "B")==0) {
	float *data = ypush_f(tscreen_handler->d_B->dims_data);
	tscreen_handler->d_B->device2host(data);
      } else if (strcmp(type_data, "istx")==0) {
	int *data = ypush_i(tscreen_handler->d_istencilx->dims_data);
	tscreen_handler->d_istencilx->device2host((unsigned int *)data);
      }  else if (strcmp(type_data, "isty")==0) {
	int *data = ypush_i(tscreen_handler->d_istencily->dims_data);
	tscreen_handler->d_istencily->device2host((unsigned int *)data);
      }  else if (strcmp(type_data, "values")==0) {
	int *data = ypush_i(tscreen_handler->d_tscreen->d_screen->dims_data);
	cutilSafeCall(cudaMemcpy((void **)&data,
				 tscreen_handler->d_tscreen->d_screen->values,
				 sizeof(unsigned int)*tscreen_handler->d_tscreen->d_screen->nb_elem,
				 cudaMemcpyDeviceToHost));
      } 
    } 
  }

  void
  Y_get_tscreen_update(int argc)
  { 

    if (yarg_subroutine())y_error("can only be called as a function");
    else {
      atmos_struct *handle = (atmos_struct *)yget_obj(argc-1,&yAtmos);
      if (handle->device != activeDevice) yoga_setDevice(handle->device);
      float alt = ygets_f(argc-2);
      yoga_atmos *atmos_handler = (yoga_atmos *)handle->yoga_atmos;
      yoga_tscreen *tscreen_handler = (yoga_tscreen *)(atmos_handler->d_screens[alt]);
      float *data = ypush_f(tscreen_handler->d_ytmp->dims_data);
      tscreen_handler->d_ytmp->device2host(data);
    } 
  }


  /*
 ____                           
/ ___|  ___  _   _ _ __ ___ ___ 
\___ \ / _ \| | | | '__/ __/ _ \
 ___) | (_) | |_| | | | (_|  __/
|____/ \___/ \__,_|_|  \___\___|
                                
   */

  typedef struct source_struct {
    void *yoga_source;
    int device;    
  } source_struct;


  void source_free(void *obj) {
    source_struct *handler = (source_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_source *source_obj_handler = (yoga_source *)(handler->yoga_source);
      delete source_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void source_print(void *obj)
  {
    source_struct *handler = (source_struct *)obj;
    yoga_source *source_obj_handler = (yoga_source *)(handler->yoga_source);
    cout << "Yoga Source Object : " << endl;
    cout << "Source @ pos. " << source_obj_handler->tposx << "\" : " <<  source_obj_handler->tposy
	 << " \"  |" << " mag : " << source_obj_handler->mag << " |"
	 << " lambda : " << source_obj_handler->lambda << " micron "<<endl;
  }


  string y_sourceName("ySource Object");
  static y_userobj_t ySource =
    {const_cast<char*>(y_sourceName.c_str()), 
     &source_free, 0, 0, 0, 0};


  /*
 _____                    _   
|_   _|_ _ _ __ __ _  ___| |_ 
  | |/ _` | '__/ _` |/ _ \ __|
  | | (_| | | | (_| |  __/ |_ 
  |_|\__,_|_|  \__, |\___|\__|
               |___/          
   */

  typedef struct target_struct {
    void *yoga_target;
    int device;    
  } target_struct;

  void target_free(void *obj) {
    target_struct *handler = (target_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_target *target_obj_handler = (yoga_target *)(handler->yoga_target);
      delete target_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void target_print(void *obj)
  {
    target_struct *handler = (target_struct *)obj;
    yoga_target *target_obj_handler = (yoga_target *)(handler->yoga_target);
    cout << "Yoga Target Object : " << endl;
    cout << "Contains " << target_obj_handler->ntargets << " target(s) : " << endl;
    for (int i=0;i<target_obj_handler->ntargets;i++) {
      cout << "Source # " << i+1 << " @ pos. " << target_obj_handler->d_targets.at(i)->tposx 
	   << "\" : " << target_obj_handler->d_targets.at(i)->tposy << " \"  |" << " mag : " 
	   << target_obj_handler->d_targets.at(i)->mag << " |" << " lambda : " 
	   << target_obj_handler->d_targets.at(i)->lambda << " micron "<<endl;
    }
  }

  string y_TargetName("yTarget Object");
  static y_userobj_t yTarget =
    {const_cast<char*>(y_TargetName.c_str()), 
     &target_free, &target_print, 0, 0, 0};

  void
  Y_yoga_target(int argc)
 {    
    long ntot;
    long dims[Y_DIMSIZE];

    try { 
      int odevice = activeDevice;
      int ntargets = ygets_i(argc-1);

      float *xpos;
      xpos = ygeta_f(argc-2, &ntot,dims);
      if (ntot != ntargets) y_error("wrong dimension for xpos");

      float *ypos;
      ypos = ygeta_f(argc-3, &ntot,dims);
      if (ntot != ntargets) y_error("wrong dimension for screens ypos");

      float *lambda;
      lambda = ygeta_f(argc-4, &ntot,dims);
      if (ntot != ntargets) y_error("wrong dimension for screens lambda");

      float *mag;
      mag = ygeta_f(argc-5, &ntot,dims);
      if (ntot != ntargets) y_error("wrong dimension for screens mag");

      long *sizes;
      sizes = ygeta_l(argc-6, &ntot,dims);

      if (argc > 6) odevice = ygets_i(argc-7);
      if (odevice != activeDevice) {
	if (odevice > _nbDevice()-1) {
	  cout << "### Invalid Device Id : " << odevice <<" Your system has only " << 
	    _nbDevice() << " device(s) available" << endl;
	  cout << "Using active Device with Id : " << activeDevice << endl;
	} else {
	  cudaDeviceProp deviceProp;
	  cutilSafeCall( cudaGetDeviceProperties(&deviceProp, odevice) );
	  yoga_setDevice(odevice);
	}
      }

      target_struct *handle=(target_struct *)ypush_obj(&yTarget, sizeof(target_struct));
      handle->device = odevice;
      
      handle->yoga_target = new yoga_target(ntargets,xpos,ypos,lambda,mag,sizes);
      
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    } catch( ... ) { 
      stringstream buf;
      buf << "unknown error with yoga_target construction in "<<__FILE__ << "@" << __LINE__ << endl;
      y_error(buf.str().c_str()); 
    }
  }
 
  void
  Y_target_addlayer(int argc)
  {    
    long ntot;
    long dims[Y_DIMSIZE];
    target_struct *handler = (target_struct *)yget_obj(argc-1,&yTarget);
    yoga_target *target_handler = (yoga_target *)(handler->yoga_target);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int ntarget = ygets_i(argc-2);

    char *type = ygets_q(argc-3);

    float alt = ygets_f(argc-4);

    float *xref = ygeta_f(argc-5, &ntot,dims);
    if (ntot != target_handler->d_targets.at(ntarget)->npos)
      y_error("wrong size for xref");
    float *yref = ygeta_f(argc-6, &ntot,dims);
    if (ntot != target_handler->d_targets.at(ntarget)->npos)
      y_error("wrong size for yref");

    target_handler->d_targets.at(ntarget)->add_layer(type,alt,xref,yref,ntot);
  }

  void
  Y_target_atmostrace(int argc)
  {    
    target_struct *handler = (target_struct *)yget_obj(argc-1,&yTarget);
    yoga_target *target_handler = (yoga_target *)(handler->yoga_target);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int ntarget = ygets_i(argc-2);

    atmos_struct *handler_a = (atmos_struct *)yget_obj(argc-3,&yAtmos);
    yoga_atmos *atmos_handler = (yoga_atmos *)(handler_a->yoga_atmos);

    target_handler->d_targets.at(ntarget)->raytrace(atmos_handler);
  }

  void
  Y_target_getphase(int argc)
  {    
    target_struct *handler = (target_struct *)yget_obj(argc-1,&yTarget);
    yoga_target *target_handler = (yoga_target *)(handler->yoga_target);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int ntarget = ygets_i(argc-2);

    float *data = ypush_f(target_handler->d_targets.at(ntarget)->d_phase->d_screen->dims_data);
    target_handler->d_targets.at(ntarget)->d_phase->d_screen->device2host(data);
  }

 void
  Y_target_getamplipup(int argc)
  {    
    target_struct *handler = (target_struct *)yget_obj(argc-1,&yTarget);
    yoga_target *target_handler = (yoga_target *)(handler->yoga_target);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int ntarget = ygets_i(argc-2);

    long * ndims_data = new long[4];
    ndims_data[0]=3;ndims_data[1]=2;
    memcpy(&ndims_data[2],&(target_handler->d_targets.at(ntarget)->d_amplipup->getDims()[1]),sizeof(long)*2);
    float *data = ypush_f(ndims_data);
    target_handler->d_targets.at(ntarget)->d_amplipup->device2host((cuFloatComplex *)data);
    delete ndims_data;
  }


  /*
 ____  _                    
|  _ \| |__   __ _ ___  ___ 
| |_) | '_ \ / _` / __|/ _ \
|  __/| | | | (_| \__ \  __/
|_|   |_| |_|\__,_|___/\___|
                            
   */

  typedef struct phase_struct {
    void *yoga_phase;
    int device;    
  } phase_struct;


  void phase_free(void *obj) {
    phase_struct *handler = (phase_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_phase *phase_obj_handler = (yoga_phase *)(handler->yoga_phase);
      delete phase_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void phase_print(void *obj)
  {
    phase_struct *handler = (phase_struct *)obj;
    yoga_phase *phase_obj_handler = (yoga_phase *)(handler->yoga_phase);
    cout << "Yoga phase Object : " << phase_obj_handler->screen_size << "x"<< 
      phase_obj_handler->screen_size << endl;
  }


  void phase_eval(void *obj, int n) {
    phase_struct *handler = (phase_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    yoga_phase *phase_handler = (yoga_phase *)(handler->yoga_phase);
    float *data = ypush_f(phase_handler->d_screen->dims_data);
    phase_handler->d_screen->device2host(data);
  }

  string y_phaseName("yPhase Object");
  static y_userobj_t yPhase =
    {const_cast<char*>(y_phaseName.c_str()), 
     &phase_free, &phase_print, &phase_eval, 0, 0};

  void
  Y_yoga_phase(int argc)
  {    
    try { 
      int odevice = activeDevice;
      long size = ygets_l(argc-1);

      if (argc > 1) odevice = ygets_i(argc-2);
      if (odevice != activeDevice) {
	if (odevice > _nbDevice()-1) {
	  cout << "### Invalid Device Id : " << odevice <<" Your system has only " << 
	    _nbDevice() << " device(s) available" << endl;
	  cout << "Using active Device with Id : " << activeDevice << endl;
	} else {
	  cudaDeviceProp deviceProp;
	  cutilSafeCall( cudaGetDeviceProperties(&deviceProp, odevice) );
	  yoga_setDevice(odevice);
	}
      }

      phase_struct *handle=(phase_struct *)ypush_obj(&yPhase, sizeof(phase_struct));
      handle->device = odevice;
      
      handle->yoga_phase = new yoga_phase(size);
      
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    } catch( ... ) { 
      stringstream buf;
      buf << "unknown error with yoga_phase construction in "<<__FILE__ << "@" << __LINE__ << endl;
      y_error(buf.str().c_str()); 
    }
  }

  void
  Y_phase_copy(int argc)
  { 

    phase_struct *handle1 = (phase_struct *)yget_obj(argc-1,&yPhase);
    phase_struct *handle2 = (phase_struct *)yget_obj(argc-2,&yPhase);

    yoga_phase *phase_handler1 = (yoga_phase *)handle1->yoga_phase;
    yoga_phase *phase_handler2 = (yoga_phase *)handle2->yoga_phase;

    phase_copy(phase_handler1,phase_handler2);

  }

  void
  Y_phase_set(int argc)
  { 
    long ntot;
    long dims[Y_DIMSIZE];
    phase_struct *handle = (phase_struct *)yget_obj(argc-1,&yPhase);
    yoga_phase *phase_handler = (yoga_phase *)handle->yoga_phase;
    float *a = ygeta_f(argc-2, &ntot,dims);
    if (ntot != phase_handler->d_screen->nb_elem)
      y_error("wrong size for array");

    phase_handler->d_screen->host2device(a);

  }


  /*
__        _______ ____  
\ \      / /  ___/ ___| 
 \ \ /\ / /| |_  \___ \ 
  \ V  V / |  _|  ___) |
   \_/\_/  |_|   |____/ 
                        
   */

  typedef struct wfs_struct {
    void *yoga_wfs;
    int device;    
  } wfs_struct;


  void wfs_free(void *obj) {
    wfs_struct *handler = (wfs_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_wfs *wfs_obj_handler = (yoga_wfs *)(handler->yoga_wfs);
      delete wfs_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void wfs_print(void *obj)
  {
    wfs_struct *handler = (wfs_struct *)obj;
    yoga_wfs *wfs_obj_handler = (yoga_wfs *)(handler->yoga_wfs);
    cout << "Yoga wfs Object : " << endl;
  }


  string y_wfsName("yWfs Object");
  static y_userobj_t yWfs =
    {const_cast<char*>(y_wfsName.c_str()), 
     &wfs_free, &wfs_print, 0, 0, 0};

  void
  Y_yoga_wfs(int argc)
  //long nxsub, long nvalid, long npix, long nrebin, long nfft
  {    
    try { 
      int odevice = activeDevice;
      long nxsub = ygets_l(argc-1);
      long nvalid = ygets_l(argc-2);
      long npix = ygets_l(argc-3);
      long nphase = ygets_l(argc-4);
      long nrebin = ygets_l(argc-5);
      long nfft = ygets_l(argc-6);
      long ntot = ygets_l(argc-7);
      long npup = ygets_l(argc-8);
      
      if (argc > 8) odevice = ygets_i(argc-9);
      if (odevice != activeDevice) {
	if (odevice > _nbDevice()-1) {
	  cout << "### Invalid Device Id : " << odevice <<" Your system has only " << 
	    _nbDevice() << " device(s) available" << endl;
	  cout << "Using active Device with Id : " << activeDevice << endl;
	} else {
	  cudaDeviceProp deviceProp;
	  cutilSafeCall( cudaGetDeviceProperties(&deviceProp, odevice) );
	  yoga_setDevice(odevice);
	}
      }
      
      wfs_struct *handle=(wfs_struct *)ypush_obj(&yWfs, sizeof(wfs_struct));
      handle->device = odevice;
      
      handle->yoga_wfs = new yoga_wfs(nxsub,nvalid,npix,nphase,nrebin,nfft,ntot,npup);
      
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    } catch( ... ) { 
      stringstream buf;
      buf << "unknown error with yoga_wfs construction in "<<__FILE__ << "@" << __LINE__ << endl;
      y_error(buf.str().c_str()); 
    }
  }

  void
  Y_wfs_initgs(int argc)
 {    
   long ntot;
   long dims[Y_DIMSIZE];
   
   wfs_struct *handle = (wfs_struct *)yget_obj(argc-1,&yWfs);
   float xpos   = ygets_f(argc-2);
   float ypos   = ygets_f(argc-3);
   float lambda = ygets_f(argc-4);
   float mag    = ygets_f(argc-5);
   long  size   = ygets_l(argc-6);
   
   if (handle->device != activeDevice) {
     yoga_setDevice(handle->device);
   }
   
   yoga_wfs *wfs_handler = (yoga_wfs *)handle->yoga_wfs;
   wfs_handler->wfs_initgs(xpos,ypos,lambda,mag,size);
 }


  typedef struct sensors_struct {
    void *yoga_sensors;
    int device;    
  } sensors_struct;


  void sensors_free(void *obj) {
    sensors_struct *handler = (sensors_struct *)obj;
    if (handler->device != activeDevice) yoga_setDevice(handler->device);
    try{
      yoga_sensors *sensors_obj_handler = (yoga_sensors *)(handler->yoga_sensors);
      delete sensors_obj_handler;
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    }
  }

  void sensors_print(void *obj)
  {
    sensors_struct *handler = (sensors_struct *)obj;
    yoga_sensors *sensors_obj_handler = (yoga_sensors *)(handler->yoga_sensors);
    cout << "Yoga sensors Object" << endl;

    cout << "Contains " << sensors_obj_handler->nsensors << " WFS(s) : " << endl;
    cout << "WFS #" << " | " << "Nsubaps" << " | " << "Nvalid" << " | " << "Npix" << " | " 
	 << "Nphase" << " | " << "Nfft" << " | " << "Nrebin" << " | " << "Npup" << endl;
 
    for (int i=0;i<sensors_obj_handler->nsensors;i++) {
      cout << "   " << i+1 << "  | " << sensors_obj_handler->d_wfs.at(i)->nxsub
	   << " x " << sensors_obj_handler->d_wfs.at(i)->nxsub << " |   "
	   << sensors_obj_handler->d_wfs.at(i)->nvalid << "   |  " 
	   << sensors_obj_handler->d_wfs.at(i)->npix << "  |   " 
	   << sensors_obj_handler->d_wfs.at(i)->nphase << "   |  "
	   << sensors_obj_handler->d_wfs.at(i)->nfft <<  "  |   "
	   << sensors_obj_handler->d_wfs.at(i)->nrebin << "   |  "
	   << sensors_obj_handler->d_wfs.at(i)->npup << endl;
    }
  }


  string y_sensorsName("ySensors Object");
  static y_userobj_t ySensors =
    {const_cast<char*>(y_sensorsName.c_str()), 
     &sensors_free, &sensors_print, 0, 0, 0};

  void
  Y_yoga_sensors(int argc)
  //long *nxsub,long *nvalid,long *npix,long *nrebin,long *nfft
 {    
    long ntot;
    long dims[Y_DIMSIZE];

    try { 
      int odevice = activeDevice;
      int nsensors = ygets_i(argc-1);

      long *nxsub = ygeta_l(argc-2, &ntot,dims);

      long *nvalid = ygeta_l(argc-3, &ntot,dims);

      long *npix = ygeta_l(argc-4, &ntot,dims);

      long *nphase = ygeta_l(argc-5, &ntot,dims);

      long *nrebin = ygeta_l(argc-6, &ntot,dims);

      long *nfft = ygeta_l(argc-7, &ntot,dims);

      long *ntota = ygeta_l(argc-8, &ntot,dims);

      long npup = ygets_l(argc-9);

      if (argc > 9) odevice = ygets_i(argc-10);
      if (odevice != activeDevice) {
	if (odevice > _nbDevice()-1) {
	  cout << "### Invalid Device Id : " << odevice <<" Your system has only " << 
	    _nbDevice() << " device(s) available" << endl;
	  cout << "Using active Device with Id : " << activeDevice << endl;
	} else {
	  cudaDeviceProp deviceProp;
	  cutilSafeCall( cudaGetDeviceProperties(&deviceProp, odevice) );
	  yoga_setDevice(odevice);
	}
      }

      sensors_struct *handle=(sensors_struct *)ypush_obj(&ySensors, sizeof(sensors_struct));
      handle->device = odevice;
      
      handle->yoga_sensors = new yoga_sensors(nsensors,nxsub,nvalid,npix,nphase,nrebin,nfft,ntota,npup);
      
    } catch ( string msg ) {
      y_error(msg.c_str());
    } catch ( char const * msg ) {
      y_error(msg);
    } catch( ... ) { 
      stringstream buf;
      buf << "unknown error with yoga_sensors construction in "<<__FILE__ << "@" << __LINE__ << endl;
      y_error(buf.str().c_str()); 
    }
  }

  void
  Y_sensors_initgs(int argc)
 {    
   long ntot;
   long dims[Y_DIMSIZE];
   
   sensors_struct *handle = (sensors_struct *)yget_obj(argc-1,&ySensors);
   yoga_sensors *sensors_handler = (yoga_sensors *)handle->yoga_sensors;
   float *xpos   = ygeta_f(argc-2, &ntot,dims);
   if (ntot != sensors_handler->nsensors) y_error("wrong dimension for xpos");
   float *ypos   = ygeta_f(argc-3, &ntot,dims);
   if (ntot != sensors_handler->nsensors) y_error("wrong dimension for ypos");
   float *lambda = ygeta_f(argc-4, &ntot,dims);
   if (ntot != sensors_handler->nsensors) y_error("wrong dimension for lambda");
   float *mag    = ygeta_f(argc-5, &ntot,dims);
   if (ntot != sensors_handler->nsensors) y_error("wrong dimension for mag");
   long  *size   = ygeta_l(argc-6, &ntot,dims);
   if (ntot != sensors_handler->nsensors) y_error("wrong dimension for size");
   
   if (handle->device != activeDevice) {
     yoga_setDevice(handle->device);
   }
   
   sensors_handler->sensors_initgs(xpos,ypos,lambda,mag,size);
 }

  void
  Y_sensors_addlayer(int argc)
  {    
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
    yoga_sensors *sensors_handler = (yoga_sensors *)(handler->yoga_sensors);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int nsensor = ygets_i(argc-2);

    char *type = ygets_q(argc-3);

    float alt = ygets_f(argc-4);

    float *xref = ygeta_f(argc-5, &ntot,dims);
    if (ntot != sensors_handler->d_wfs.at(nsensor)->d_gs->npos)
      y_error("wrong size for xref");
    float *yref = ygeta_f(argc-6, &ntot,dims);
    if (ntot != sensors_handler->d_wfs.at(nsensor)->d_gs->npos)
      y_error("wrong size for yref");

    sensors_handler->d_wfs.at(nsensor)->d_gs->add_layer(type,alt,xref,yref,ntot);
  }

  void
  Y_sensors_initarr(int argc)
  {    
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
    yoga_sensors *sensors_handler = (yoga_sensors *)(handler->yoga_sensors);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int nsensor = ygets_i(argc-2);

    int *phasemap = ygeta_i(argc-3, &ntot,dims);
    int *hrmap = ygeta_i(argc-4, &ntot,dims);
    int *imamap = ygeta_i(argc-5, &ntot,dims);
    int *binmap = ygeta_i(argc-6, &ntot,dims);
    float *offsets = ygeta_f(argc-7, &ntot,dims);
    float *pupil = ygeta_f(argc-8, &ntot,dims);

    sensors_handler->d_wfs.at(nsensor)->wfs_initarrays(phasemap,hrmap,imamap,binmap,offsets,pupil);
  }

  void
  Y_sensors_compimg(int argc)
  {    
    long ntot;
    long dims[Y_DIMSIZE];
    sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
    yoga_sensors *sensors_handler = (yoga_sensors *)(handler->yoga_sensors);

    int nsensor = ygets_i(argc-2);

    atmos_struct *handlera = (atmos_struct *)yget_obj(argc-3,&yAtmos);
    yoga_atmos *atmos_handler = (yoga_atmos *)(handlera->yoga_atmos);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    sensors_handler->d_wfs.at(nsensor)->comp_image(atmos_handler,handler->device);
  }

  void
  Y_sensors_getimg(int argc)
  {    
    sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
    yoga_sensors *sensors_handler = (yoga_sensors *)(handler->yoga_sensors);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int nsensor = ygets_i(argc-2);

    float *data = ypush_f(sensors_handler->d_wfs.at(nsensor)->d_binimg->dims_data);
    sensors_handler->d_wfs.at(nsensor)->d_binimg->device2host(data);
  }

  void
  Y_sensors_getdata(int argc)
  {    
    sensors_struct *handler = (sensors_struct *)yget_obj(argc-1,&ySensors);
    yoga_sensors *sensors_handler = (yoga_sensors *)(handler->yoga_sensors);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int nsensor = ygets_i(argc-2);

    char *type_data = ygets_q(argc-3);
    if (strcmp(type_data, "amplipup")==0) {
      long *ndims_data = new long[5];
      ndims_data[0]=4;ndims_data[1]=2;
      memcpy(&ndims_data[2],&(sensors_handler->d_wfs.at(nsensor)->d_camplipup->getDims()[1]),sizeof(long)*3);
      float *data = ypush_f(ndims_data);
      sensors_handler->d_wfs.at(nsensor)->d_camplipup->device2host((cuFloatComplex *)data);
    }
    if (strcmp(type_data, "amplifoc")==0) {
      long *ndims_data = new long[5];
      ndims_data[0]=4;ndims_data[1]=2;
      memcpy(&ndims_data[2],&(sensors_handler->d_wfs.at(nsensor)->d_camplifoc->getDims()[1]),sizeof(long)*3);
      float *data = ypush_f(ndims_data);
      sensors_handler->d_wfs.at(nsensor)->d_camplifoc->device2host((cuFloatComplex *)data);
    }
    if (strcmp(type_data, "hrimg")==0) {
      float *data = ypush_f(sensors_handler->d_wfs.at(nsensor)->d_hrimg->dims_data);
      sensors_handler->d_wfs.at(nsensor)->d_hrimg->device2host(data);
    }
    if (strcmp(type_data, "bincube")==0) {
      float *data = ypush_f(sensors_handler->d_wfs.at(nsensor)->d_bincube->dims_data);
      sensors_handler->d_wfs.at(nsensor)->d_bincube->device2host(data);
    }
    if (strcmp(type_data, "phase")==0) {
      float *data = ypush_f(sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->dims_data);
      sensors_handler->d_wfs.at(nsensor)->d_gs->d_phase->d_screen->device2host(data);
    }
    if (strcmp(type_data, "totimg")==0) {
      float *data = ypush_f(sensors_handler->d_wfs.at(nsensor)->d_totimg->dims_data);
      sensors_handler->d_wfs.at(nsensor)->d_totimg->device2host(data);
    }

  }

  /*
  ____ _       _           _ 
 / ___| | ___ | |__   __ _| |
| |  _| |/ _ \| '_ \ / _` | |
| |_| | | (_) | |_) | (_| | |
 \____|_|\___/|_.__/ \__,_|_|
                             
   */

  void
  Y_move_atmos(int argc)
  {    
    atmos_struct *handler_a = (atmos_struct *)yget_obj(argc-1,&yAtmos);
    yoga_atmos *atmos_handler = (yoga_atmos *)(handler_a->yoga_atmos);

    if (handler_a->device != activeDevice) yoga_setDevice(handler_a->device);

    move_atmos(atmos_handler);
  }

  void
  Y_move_sky(int argc)
  {    
    atmos_struct *handler_a = (atmos_struct *)yget_obj(argc-1,&yAtmos);
    yoga_atmos *atmos_handler = (yoga_atmos *)(handler_a->yoga_atmos);

    target_struct *handler = (target_struct *)yget_obj(argc-2,&yTarget);
    yoga_target *target_handler = (yoga_target *)(handler->yoga_target);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    move_atmos(atmos_handler,target_handler);
  }

  void
  Y_target_getimage(int argc)
  {    
    target_struct *handler = (target_struct *)yget_obj(argc-1,&yTarget);
    yoga_target *target_handler = (yoga_target *)(handler->yoga_target);

    atmos_struct *handler_a = (atmos_struct *)yget_obj(argc-2,&yAtmos);
    yoga_atmos *atmos_handler = (yoga_atmos *)(handler_a->yoga_atmos);

    if (handler->device != activeDevice) yoga_setDevice(handler->device);

    int ntarget = ygets_i(argc-3);

    float *data = ypush_f(target_handler->d_targets.at(ntarget)->d_image->dims_data);

    target_handler->d_targets.at(ntarget)->raytrace(atmos_handler);
    target_handler->d_targets.at(ntarget)->comp_image(atmos_handler->d_pupil->d_data);
    target_handler->d_targets.at(ntarget)->d_image->device2host(data);
  }

}

