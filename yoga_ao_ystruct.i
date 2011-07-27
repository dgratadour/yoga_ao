struct geom_struct
{
  long  ssize;       // linear size of full image (in pixels)
  float zenithangle; // observations zenith angle (in deg)
  
  // internal keywords
  long  pupdiam;     // linear size of total pupil (in pixels)
  float cent;        // central point of the simulation
  pointer _ipupil;   // total pupil (include full guard band)
  pointer _mpupil;   // medium pupil (part of the guard band)
  pointer _spupil;   // small pupil (without guard band)
  long  _p1;         // min x,y for valid points in mpupil
  long  _p2;         // max x,y for valid points in mpupil
  long  _n;          // linear size of mpupil
  long  _n1;         // max x,y for valid points in ipupil
  long  _n2;         // min x,y for valid points in ipupil
};

struct tel_struct
{
  float diam;        // telescope diameter (in meters)
  float cobs;        // central obstruction ratio
};

struct atmos_struct
{
  long    nscreens;    // number of turbulent layers
  float   r0;          // global r0 @ 0.5µm
  float   pupixsize;   // pupil piwel size (in meters)
  pointer dim_screens; // linear size of phase screens
  pointer alt;         // altitudes of each layer
  pointer winddir;     // wind directions of each layer
  pointer windspeed;   // wind speeds of each layer
  pointer frac;        // fraction of r0 for each layer
  pointer deltax;      // x translation speed (in pix / iteration) for each layer
  pointer deltay;      // y translation speed (in pix / iteration) for each layer
};

struct target_struct
{
  long    ntargets;  // number of targets
  pointer lambda;    // observation wavelength for each target
  pointer xpos;      // x positions on sky (in arcsec) for each target
  pointer ypos;      // y positions on sky (in arcsec) for each target
  pointer mag;       // magnitude for each target
};

struct wfs_struct
{
  long  nxsub;          // linear number of subaps
  long  npix;           // number of pixels per subap
  float pixsize;        // pixel size (in arcsec) for a subap
  float lambda;         // observation wavelength (in µm) for a subap
  float optthroughput;  // wfs global throughput
  float fracsub;        // minimal illumination fraction for valid subaps
  
  //target kwrd
  float xpos;      // guide star x position on sky (in arcsec) 
  float ypos;      // guide star x position on sky (in arcsec) 
  float gsalt;     // altitude of guide star (in m) 0 if ngs 
  float gsmag;     // magnitude of guide star
  float zerop;     // detector zero point
  
  // lgs only
  float lgsreturnperwatt;  // return per watt factor (high season : 10 ph/cm2/s/W)
  float laserpower;        // laser power in W
  float lltx;              // x position (in meters) of llt
  float llty;              // y position (in meters) of llt
  string proftype;         // type of sodium profile "gauss", "exp", etc ...
  float beamsize;          // laser beam fwhm on-sky (in arcsec)

  //internal kwrd
  long  _pdiam;          // pupil diam for a subap (in pixel)
  long  _Nfft;           // array size for fft for a subap (in pixel)
  long  _Ntot;           // total size of hr image for a subap (in pixel)
  long  _nrebin;         // rebin factor from hr to binned image for a subap 
  long  _nvalid;         // number of valid subaps

  float _nphotons;       // number of photons per subap
  float _qpixsize;       // quantum pixel size for the simulation
  
  pointer _istart;       // x start indexes for cutting phase screens 
  pointer _jstart;       // y start indexes for cutting phase screens 
  pointer _validsubs;    // (i,j) indices of valid subaps
  pointer _isvalid;      // array of 0/1 for valid subaps
  pointer _phasemap;     // array of pixels transform from phase screen into
                         // subaps phase screens
  pointer _hrmap;        // array of pixels transform from minimal FoV image to
                         // full FoV image (array of 0 if the same)
  pointer _binmap;       // array of pixels transform from full FoV hr images to
                         // binned images
  pointer _imamap;       // array of pixels transform from subaps binned images to
                         // total wfs image
  pointer _halfxy;       // phase offset for 1/2 pixel shift in (x,y)

  pointer xprofbin;      // pixel transform from hr profile to binned profile
  pointer deltaprofbin;  // pixel transform from hr profile to binned profile
  pointer _lgskern;      // lgs kernels for each subap
};

struct loop_struct
{
  long  niter;      // number of iterations
  float ittime;     // iteration time (in sec)
};

