#ifndef _YOGA_WFS_H_
#define _YOGA_WFS_H_

#include <vector>
#include <map>
#include <yoga_target.h>
#include <yoga_phase.h>

using namespace std;

class yoga_wfs {
 public:
  long                     nxsub;
  long                     nvalid;
  long                     npix;
  long                     nrebin;
  long                     nfft;
  long                     ntot;
  long                     npup;
  long                     nphase;
  bool                     lgs;

  yoga_source              *d_gs;

  yoga_obj<cuFloatComplex> *d_camplipup;
  yoga_obj<cuFloatComplex> *d_camplifoc;
  yoga_obj<float>          *d_hrimg;
  yoga_obj<float>          *d_totimg;
  yoga_obj<float>          *d_bincube;
  yoga_obj<float>          *d_binimg;
  yoga_obj<int>            *d_phasemap;
  yoga_obj<int>            *d_hrmap;
  yoga_obj<int>            *d_binmap;
  yoga_obj<int>            *d_imamap;
  yoga_obj<float>          *d_offsets;
  yoga_obj<float>          *d_pupil;
  yoga_obj<float>          *d_lgskern;
  yoga_obj<cuFloatComplex> *d_ftlgskern;
  yoga_obj<cuFloatComplex> *d_fttotim;

  yoga_obj<float>          *d_slopes;     
  yoga_obj<float>          *d_refslopes;

 public:
  yoga_wfs(long nxsub, long nvalid, long npix, long nphase, long nrebin, long nfft, long ntot, long npup,int lgs);
  yoga_wfs(const yoga_wfs& wfs);
  ~yoga_wfs();

  int wfs_initarrays(int *phasemap,int *hrmap,int *imamap,int *binmap,float *offsets, float *pupil);
  int wfs_initgs(float xpos,float ypos,float lambda, float mag, long size);
  int comp_image(yoga_atmos *yatmos, int device);
  int load_kernels(float *lgskern, int device);
};

class yoga_sensors {
 public:
  int                 nsensors;
  vector<yoga_wfs *>  d_wfs;
     
 public:
  yoga_sensors(int nwfs,long *nxsub,long *nvalid,long *npix,long *nphase, long *nrebin,long *nfft, long *ntot, int *lgs, long npup);
  ~yoga_sensors();

  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size);
};

int fillcamplipup(cuFloatComplex *amplipup, float *phase,float *offset, float *mask, int *indx, int Nfft, int Npup, int Nsub, int npup, int device);
int fillbincube(float *bcube, float *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin, int Nsub, int device);
int fillbinimg(float *bimage, float *bcube, int *indxbin, int Npix, int Nsub, int device);
int indexfill(float *d_odata,float *d_idata,int *indx,int nx, int Nx,int N,int device);
int launch_cfillrealp(cuFloatComplex *d_odata,float *d_idata,int N,int device);
int launch_conv_krnl(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device);
int launch_cgetrealp(float *d_odata,cuFloatComplex *d_idata,int N,int device);

#endif // _YOGA_WFS_H_

