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
/* DOCUMENT yoga_atmos
   obj = yoga_atmos(nscreens,r0,size,size2,alt,wspeed,wdir,deltax,deltay,pupil[,ndevice])

   creates an yAtmos object on the gpu
   nscreens : # of layers
   r0       : r0 for each layer
   size     : linear size of the screen
   size2    : second dim of extrude matrix A : size x size2
   alt      : array of altitudes per layers
   wspeed   : array of wind speeds per layers
   wdir     : array of wind directions per layers
   deltax   : array of x displacement per iteration (one per layers)
   deltay   : array of y displacement per iteration (one per layers)
   pupil    : array containing the pupil
   
   SEE ALSO:
 */
extern init_tscreen;
/* DOCUMENT init_tscreen
   init_tscreen,yoga_atmos_obj,altitude,a,b,istencilx,istencily,seed
     
   loads on the gpu in an yAtmos object and for a given screen data needed for extrude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   a              : matrix A
   b              : matrix B
   istencilx      : stencil for x direction
   istencily      : stencil for y direction
   seed           : seed for random numbers
   
   SEE ALSO:
 */
extern get_tscreen;
/* DOCUMENT get_tscreen
   screen = get_tscreen(yoga_atmos_obj,altitude)
     
   returns the screen in an yAtmos object and for a given altitude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   
   SEE ALSO:
 */
extern get_spupil;
/* DOCUMENT get_spupil
   pup = get_spupil(yoga_atmos_obj)
     
   returns the pupil in an yAtmos object
   yoga_atmos_obj : the yAtmos object
   
   SEE ALSO:
 */
extern get_tscreen_config;
/* DOCUMENT get_tscreen_config
   arr = get_tscreen_config(yoga_atmos_obj,altitude,str)
     
   returns config data in an yAtmos object and for a given altitude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   str            : type of data "A" or "B" or "istx" or "isty" or "values"
   
   SEE ALSO:
 */
extern get_tscreen_update;
/* DOCUMENT get_tscreen_update
   vect = get_tscreen_update(yoga_atmos_obj,altitude)
     
   returns only the update vector in an yAtmos object and for a given altitude
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   
   SEE ALSO:
 */
extern extrude_tscreen;
/* DOCUMENT extrude_tscreen
   extrude_tscreen,yoga_atmos_obj,altitude[,dir]
     
   executes one col / row screen extrusion for a given altitude in an yAtmos object 
   yoga_atmos_obj : the yAtmos object
   altitude       : the altitude of the given screen
   dir            : 0 (default) = column / 1 = row
   
   SEE ALSO:
 */

// targets
extern yoga_target;
/* DOCUMENT yoga_target
   obj = yoga_target(ntargets,xpos,ypos,lambda,mag,sizes[,ndevice])
     
   creates an yTarget object on the gpu
   ntargets : # of targets
   xpos     : array of targets x positions
   ypos     : array of targets y positions
   lambda   : array of wavelengths
   mag      : array of magnitudes
   sizes    : array of linear # of pixels in the pupil
   
   SEE ALSO:
 */
extern target_addlayer;
/* DOCUMENT target_addlayer
   target_addlayer,yoga_target_obj,ntarget,type,alt,xref,yref
   with type = "wfs" or "img"
     
   adds a layer of disturbances for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target
   type            : type of layer "atmos" or "optics"
   alt             : altitude of layer
   xref            : array of x positions of light rays in the pupil
   yref            : array of y positions of light rays in the pupil
   
   SEE ALSO:
 */
extern target_atmostrace;
/* DOCUMENT target_atmostrace
   target_atmostrace,yoga_target_obj,ntarget
     
   does raytracing through the atmos layers for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target
   
   SEE ALSO:
 */
extern target_getphase;
/* DOCUMENT target_getphase
   screen = target_getphase(yoga_target_obj,ntarget)
     
   returns the phase for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */
extern target_getamplipup;
/* DOCUMENT target_getamplipup
   ampli = target_getamplipup(yoga_target_obj,ntarget)
     
   returns the complex amplitude in the pupil plane for a given target in an yTarget object
   yoga_target_obj : the yTarget object
   ntarget         : index of the given target

   SEE ALSO:
 */

// phase
extern yoga_phase;
extern phase_copy;
extern phase_set;

// wfs
extern yoga_wfs;
/* DOCUMENT yoga_wfs
   obj = yoga_wfs(nxsub,nvalid,npix,nphase,nrebin,nfft,ntot,lgs,npup[,ndevice])
     
   creates an yWfs object on the gpu
   nxsub  : linear # of subaps
   nvalid : number of valid subaps
   npix   : linear number of cam pixels in a subap
   nphase : linear size of phase screens in a subap
   nrebin : rebin factor from hr img to cam pixels
   nfft   : linear size of fft arrays for psf computations for a subap
   ntot   : linear size of hr image (total FoV) for a subap
   lgs    : flag for lgs mode : 1 if lgs wfs
   npup   : linear size of total pupil image
      
   SEE ALSO:
 */
extern wfs_initgs;
/* DOCUMENT wfs_initgs
   wfs_initgs,yoga_wfs_obj,xpos,ypos,lambda,mag,size
     
   inits the guide star for an yWfs object
   yoga_wfs_obj : the yWfs object
   xpos         : x position of the guide star
   ypos         : y position of the guide star
   lambda       : observing wavelength for the guide star
   mag          : magnitude of the guide star
   size         : linear size of total image
   
   SEE ALSO:
 */
extern yoga_sensors;
/* DOCUMENT yoga_sensors
   obj = yoga_sensors(nsensors,nxsub,nvalid,npix,nphase,nrebin,nfft,ntot,lgs,npup[,ndevice])
     
   creates an ySensors object on the gpu
   nsensors : # of wfs
   nxsub    : array of linear # of subaps
   nvalid   : array of numbers of valid subaps
   npix     : array of linear numbers of cam pixels in a subap
   nphase   : array of linear sizes of phase screens in a subap
   nrebin   : array of rebin factors from hr img to cam pixels
   nfft     : array of linear sizes of fft arrays for psf computations for a subap
   ntot     : array of linear sizes of hr image (total FoV) for a subap
   lgs      : array of flags for lgs mode : 1 if lgs wfs
   npup     : array of linear sizes of total pupil image
      
   SEE ALSO:
 */
extern sensors_initgs;
/* DOCUMENT sensors_initgs
   sensors_initgs,yoga_sensors_obj,xpos,ypos,lambda,mag,size
     
   inits the guide stars for an ySensors object
   yoga_sensors_obj : the ySensors object
   xpos             : array of x positions of the guide stars
   ypos             : array of y positions of the guide stars
   lambda           : array of observing wavelengths for the guide stars
   mag              : array of magnitudes of the guide stars
   size             : array of linear sizes of total images
   
   SEE ALSO:
 */
extern sensors_addlayer;
/* DOCUMENT sensors_addlayer
   sensors_addlayer,yoga_sensors_obj,nsensor,type,alt,xref,yref
   type not used can be "atmos"
     
   add a disturbance layer for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   type             : type of layer "atmos" or "optics"
   alt              : altitude of layer
   xref             : array of x positions of light rays in the pupil
   yref             : array of y positions of light rays in the pupil
  
   SEE ALSO:
 */
extern sensors_initarr;
/* DOCUMENT sensors_initarr
   sensors_initarr,yoga_sensors_obj,nsensor,phasemap,hrmap,imamap,binmap,offsets,pupil
     
   init arrays for image computation for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   phasemap         : array of pixels transform from phase screen into
                      subaps phase screens
   hrmap            : array of pixels transform from minimal FoV image to
                      full FoV image (array of 0 if the same)
   imamap           : array of pixels transform from subaps binned images to
                      total wfs image
   binmap           : array of pixels transform from full FoV hr images to
                      binned images
   offsets          : array of pixels offsets for subaps phase screens
   pupil            : the pupil array
   
   SEE ALSO:
 */
extern sensors_compimg;
/* DOCUMENT sensors_compimg
   sensors_compimg,yoga_sensors_obj,nsensor,yoga_atmos_obj
     
   image computation for a given sensor in a ySensors object and given turbulence in an yAtmos object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   yoga_atmos_obj   : the yAtmos object

   SEE ALSO:
 */
extern sensors_getimg;
/* DOCUMENT sensors_getimg
   img = sensors_getimg(yoga_sensors_obj,nsensor)
     
   returns the total image for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs

   SEE ALSO:
 */
extern sensors_getdata;
/* DOCUMENT sensors_getdata
   arr = sensors_getdata(yoga_sensors_obj,nsensor,type)
   type = 
     
   returns the configuration data for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   type             : type of data "amplipup" or "amplifoc" or "hrimg" or
                      "bincube" or "phase" or "totimg" or "lgskern" or "ftlgskern"

   SEE ALSO:
 */
extern sensors_loadkernels;
/* DOCUMENT sensors_loadkernels
   sensors_loadkernels,yoga_sensors_obj,nsensor,kernels
     
   load lgs kernels for a given sensor in a ySensors object
   yoga_sensors_obj : the ySensors object
   nsensor          : index of given wfs
   kernels          : array of lgs kernels for each subap

   SEE ALSO:
 */

// global
extern move_atmos;
/* DOCUMENT move_atmos
   move_atmos,yoga_atmos_obj
     
   multi-layer extrude process for a yAtmos object

   SEE ALSO:
 */
extern move_sky;
/* DOCUMENT move_sky
   move_sky,yoga_atmos_obj,yoga_target_obj
     
   multi-layer & multi-target extrude process for a yAtmos object and a yTarget object

   SEE ALSO:
 */
extern target_getimage;
/* DOCUMENT target_getimage
   img = target_getimage(yoga_target_obj,yoga_atmos_obj,ntarget)
     
   get given target image for a yTarget object and a yAtmos object
   yoga_target_obj : the yTarget object
   yoga_atmos_obj  : the yAtmos object
   ntarget         : index of given target

   SEE ALSO:
 */
