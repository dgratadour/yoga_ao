struct geom_struct
{
  long  pupdiam;
  long  ssize;
  float zenithangle;
  float cent;
  
  // internal keywords
  pointer _ipupil;
  pointer _mpupil;
  pointer _spupil;
  long  _p1;
  long  _p2;
  long  _n;
  long  _n1;
  long  _n2;
};

struct tel_struct
{
  float diam;
  float cobs;
};

struct atmos_struct
{
  long    nscreens;
  float   r0;
  float   pupixsize;
  pointer dim_screens;
  pointer alt;
  pointer winddir;
  pointer windspeed;
  pointer frac;
  pointer deltax;
  pointer deltay;
};

struct target_struct
{
  long    ntargets;
  pointer lambda;
  pointer xpos;
  pointer ypos;
  pointer mag;
};

struct wfs_struct
{
  long  nxsub;
  long  npix;
  float pixsize;
  float lambda;
  float optthroughput;
  float fracsub;
  
  //target kwrd
  float xpos;
  float ypos;
  float gsalt;
  float gsmag;
  float zerop;

  // lgs only
  float lgsreturnperwatt;
  float laserpower;
  float lltx;
  float llty;
  string proftype;

  //internal kwrd
  long  _pdiam;
  long  _Nfft;
  long  _Ntot;
  long  _nrebin;
  long  _nvalid;

  float _nphotons;
  float _qpixsize;
  
  pointer _istart;
  pointer _jstart;
  pointer _validsubs;
  pointer _isvalid;
  pointer _phasemap;
  pointer _hrmap;
  pointer _binmap;
  pointer _imamap;
  pointer _halfxy;
  pointer _subx;
  pointer _suby;

  pointer xprofbin
  pointer deltaprofbin
};

struct loop_struct
{
  long  niter;
  float ittime;
};

