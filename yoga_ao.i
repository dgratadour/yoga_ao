// Environment check
yoga_ao_top = get_env("YOGA_AO_TOP");
if (!yoga_ao_top) error,"YOGA_AO_TOP is not defined!";

mypath = anyof(split_path(get_path())==(yoga_ao_top+"/")) ? [] : get_path()+":"+yoga_ao_top+"/";
if (mypath != []) set_path,mypath;

YOGA_AO_SAVEPATH = yoga_ao_top+"/data/";
mkdirp,YOGA_AO_SAVEPATH;
// creates data dir if does not exists (else mkdirp does nothing)

require,"yoga_aolib.i"
require,"yoga_ao_ystruct.i"
require,"yoga_ao_utils.i"
require,"yoga_turbu.i"
require,"yoga_wfs.i"

func read_parfile(filename)
/* DOCUMENT read_parfile
   read_parfile,filename
     
   reads yoga parameter file filename

   SEE ALSO:
 */
{
  extern y_geom,y_atmos,y_tel,y_target,y_loop,y_wfs;

  if (!fileExist(filename)) {
    exit,swrite(format="Could not find parameter file %s !\n",filename);}

  y_geom   = geom_struct();
  y_atmos  = atmos_struct();
  y_tel    = tel_struct();
  y_target = target_struct();
  y_loop   = loop_struct();
  
  require,filename;
}


func geom_init(pupdiam)
/* DOCUMENT geom_init
   geom_init,pupdiam
     
   inits simulation geometry, depending on pupdiam
   the linear number of pixels in the pupil
   
   requires 3 externals : y_atmos, y_geom and g_atmos
   y_atmos : a y_struct for the atmosphere
   y_geom  : a y_struct for the geometry of the simulation
   g_atmos : a yAtmos object on the gpu

  SEE ALSO:
 */
{
  extern y_atmos,y_geom,g_atmos;

  if (y_geom == []) y_geom   = geom_struct();
  y_geom.pupdiam = pupdiam;
  
  y_geom.ssize  = long(2^ceil(log(y_geom.pupdiam)/log(2)+1));

  // FIX ME ! using images centered on 1/2 pixels
  y_geom.cent   = y_geom.ssize / 2 + 0.5;
  
  _p        = y_geom.pupdiam;
  _p1       = long(ceil(y_geom.cent-_p/2.));
  _p2       = long(floor(y_geom.cent+_p/2.));
  _p        = _p2-_p1+1;

  y_geom.pupdiam = _p;
  y_geom._p1     = _p1;
  y_geom._p2     = _p2;

  y_geom._n      = _p+4;
  y_geom._n1     = _p1-2;
  y_geom._n2     = _p2+2;

  // Initialize pupil array:

  ipupil = float(make_pupil(y_geom.ssize,y_geom.pupdiam,xc=y_geom.cent,yc=y_geom.cent,\
                            cobs=y_tel.cobs));

  spupil = float(make_pupil(y_geom.pupdiam,y_geom.pupdiam,xc=y_geom.pupdiam/2+0.5,yc=y_geom.pupdiam/2+0.5,\
                            cobs=y_tel.cobs));

  mpupil = float(make_pupil(y_geom._n,y_geom.pupdiam,xc=y_geom._n/2+0.5,yc=y_geom._n/2+0.5,\
                            cobs=y_tel.cobs));

  y_geom._ipupil = &ipupil;
  y_geom._spupil = &spupil;
  y_geom._mpupil = &mpupil;
}

func atmos_init(void)
/* DOCUMENT atmos_init
   atmos_init
     
   inits a yAtmos object on the gpu
   no input parameters
   
   requires 2 externals + 2 optional : y_atmos and y_geom + y_target and y_wfs
   y_atmos  : a y_struct for the atmosphere
   y_geom   : a y_struct for the geometry of the simulation
   y_target : a y_struct for the targets
   y_wfs    : a y_struct for the sensors
   creates 1 external :
   g_atmos : a yAtmos object on the gpu

  SEE ALSO:
 */
{
  extern g_atmos;
  
  y_atmos.alt = &(*y_atmos.alt / cos(y_geom.zenithangle*dtor));
  y_atmos.pupixsize = y_tel.diam / y_geom.pupdiam; // pixel size in meter

  if ((y_wfs != []) && (y_target != []))
    max_size = max(_(abs(*y_target.xpos,*y_target.ypos),abs(y_wfs.xpos,y_wfs.ypos)));
  else {
    if (y_target != [])
      max_size = max(abs(*y_target.xpos,*y_target.ypos));
    else if (y_wfs != [])
      max_size = max(abs(y_wfs.xpos,y_wfs.ypos));
    else max_size = 0.;
  }
  
  patch_diam = long(y_geom._n+2.5*max_size*4.848e-6*(*y_atmos.alt)/y_atmos.pupixsize);
  patch_diam += patch_diam % 2;

  y_atmos.dim_screens = &patch_diam;

  deltax  = y_geom.pupdiam/y_tel.diam*(*y_atmos.windspeed)*y_loop.ittime;
  deltay  = deltax*sin(dtor*(*y_atmos.winddir));
  deltax  = deltax*cos(dtor*(*y_atmos.winddir));

  y_atmos.deltax = &float(deltax);
  y_atmos.deltay = &float(deltay);

  g_atmos = yoga_atmos_create(y_atmos.nscreens,y_atmos.r0,y_atmos.pupixsize,*y_atmos.dim_screens,
                            *y_atmos.frac,*y_atmos.alt,*y_atmos.windspeed,*y_atmos.winddir,
                            *y_atmos.deltax,*y_atmos.deltay,*y_geom._spupil);

}

func wfs_init(void)
/* DOCUMENT wfs_init
   wfs_init
     
   inits a ySeznsors object on the gpu
   no input parameters
   
   requires 2 externals : y_geom +and y_wfs
   y_geom   : a y_struct for the geometry of the simulation
   y_wfs    : a y_struct for the sensors
   creates 1 external :
   g_wfs    : a ySensors object on the gpu

  SEE ALSO:
 */
{
  extern y_geom;
  extern g_wfs;

  // first get the wfs with max # of subaps
  // we'll derive the geometry from the requirements in terms of sampling
  indmax = wheremax(y_wfs.nxsub)(1);
  init_wfs_geom,indmax;
  
  // do the same for other wfs
  for (i=1;i<=numberof(y_wfs);i++) {
    if (i != wheremax(y_wfs.nxsub)(1)) {
      init_wfs_geom,i;
    }
  }
  g_wfs = yoga_sensors(numberof(y_wfs),y_wfs.nxsub,y_wfs._nvalid,y_wfs.npix,y_wfs._pdiam,
                       y_wfs._nrebin,y_wfs._Nfft,y_wfs._Ntot,y_geom._n,y_wfs._subapd,
                       y_wfs._nphotons,y_wfs.gsalt > 0);
  
  sensors_initgs,g_wfs,y_wfs.xpos,y_wfs.ypos,y_wfs.lambda,y_wfs.gsmag,
    (y_geom._n)(-::numberof(y_wfs)-1),y_wfs.noise;
  
  for (i=1;i<=numberof(y_wfs);i++) {
    sensors_initarr,g_wfs,i-1,int(*y_wfs(i)._phasemap),int(*y_wfs(i)._hrmap),
      int(*y_wfs(i)._binmap),float(*y_wfs(i)._halfxy),float(*y_geom._mpupil),
      (*y_wfs(i)._fluxPerSub)(where(*y_wfs(i)._isvalid)),int(*y_wfs(i)._isvalid),
      int((*y_wfs(i)._validsubs)(1,)-1),int((*y_wfs(i)._validsubs)(2,)-1),int(*y_wfs(i)._istart+1),int(*y_wfs(i)._jstart+1);
  }

  // lgs case
  for (cc=1;cc<=numberof(y_wfs);cc++) {
    if (y_wfs(cc).gsalt > 0) {
      // lgs mode requested
      if (y_wfs(cc).proftype == "None") y_wfs(cc).proftype = "Gauss1";
      if (y_wfs(cc).proftype == "Gauss1") profilename = "allProfileNa_withAltitude_1Gaussian.fits";
      if (y_wfs(cc).proftype == "Gauss2") profilename = "allProfileNa_withAltitude_2Gaussians.fits";
      if (y_wfs(cc).proftype == "Gauss3") profilename = "allProfileNa_withAltitude_3Gaussians.fits";
      if (y_wfs(cc).proftype == "Exp") profilename = "allProfileNa_withAltitude.fits";
      prof=fits_read(YOGA_AO_SAVEPATH+profilename);
      h    = prof(,1);
      prof = prof(,2:)(,avg);
      prep_lgs_prof,cc,prof,h,y_wfs(cc).beamsize;
      //sensors_loadkernels,g_wfs,cc-1,float(*y_wfs(cc)._lgskern);
    }
  }
}

func target_init(void)
/* DOCUMENT target_init
   target_init
     
   inits a yTarget object on the gpu
   no input parameters
   
   requires 4 externals + 1 optional : y_geom y_atmos and y_target + y_wfs
   y_geom   : a y_struct for the geometry of the simulation
   y_atmos  : a y_struct for the atmosphere
   y_target : a y_struct for the targets
   y_wfs    : a y_struct for the sensors
   creates 1 external :
   g_target    : a yTarget object on the gpu

  SEE ALSO:
 */
{
  extern g_target;
  
  type = "atmos";

  if (y_target != []) {
    sizes = y_geom.pupdiam;
    sizes = sizes(-::y_target.ntargets-1);
    
    g_target = yoga_target(y_target.ntargets,*y_target.xpos,*y_target.ypos,*y_target.lambda,*y_target.mag,sizes);

    for (cc=1;cc<=y_target.ntargets;cc++) {
      for (dd=1;dd<=y_atmos.nscreens;dd++) {
        xoff = (*y_target.xpos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
        yoff = (*y_target.ypos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
        xoff = float(xoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
        yoff = float(yoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
        target_addlayer,g_target,cc-1,type,(*y_atmos.alt)(dd),xoff,yoff;
      }
    }
  }
  
  if (y_wfs != []) {
    if ((y_wfs != []) && (g_wfs != [])) {
      for (cc=1;cc<=numberof(y_wfs);cc++) {
        for (dd=1;dd<=y_atmos.nscreens;dd++) {
          xoff = (y_wfs.xpos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
          yoff = (y_wfs.ypos)(cc)*4.848e-6*(*y_atmos.alt)(dd)/y_atmos.pupixsize;
          xoff = float(xoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
          yoff = float(yoff+((*y_atmos.dim_screens)(dd)-y_geom._n)/2);
          sensors_addlayer,g_wfs,cc-1,type,(*y_atmos.alt)(dd),xoff,yoff;
        }
      }
    }
  }
}


/*
read_parfile,"yoga_ao_test.par";
"WFS inits";
wfs_init;
"Atmos inits";
atmos_init;
"Targets inits";
target_init;
*/

/*
// check turbulence
AB, y_geom.pupdiam, A, B, ist;   // initialisation for A and B matrices for phase extrusion
p=array(0.0,y_geom.pupdiam,y_geom.pupdiam);   // init of first phase screen
for(i=1;i<=4*y_geom.pupdiam;i++) p=extrude(p, y_atmos.r0/y_atmos.pupixsize, A, B, ist);  // iterates extrusion in order
for(i=1;i<=2*y_geom.pupdiam;i++) move_atmos,g_atmos;
target_atmostrace,g_target,0,g_atmos;
res = target_getphase(g_target,0);
*/

/*
sensors_compimg,g_wfs,0,g_atmos;
res=sensors_getimg(g_wfs,0);
pli,res;
*/
