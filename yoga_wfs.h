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
  float                    subapd;
  float                    nphot;
  float                    noise;
  bool                     lgs;

  yoga_source              *d_gs;

  yoga_obj<cuFloatComplex> *d_camplipup;
  yoga_obj<cuFloatComplex> *d_camplifoc;
  yoga_obj<float>          *d_hrimg;
  yoga_obj<float>          *d_totimg;
  yoga_obj<float>          *d_bincube;
  yoga_obj<float>          *d_binimg;
  yoga_obj<float>          *d_subsum;
  yoga_obj<int>            *d_phasemap;
  yoga_obj<int>            *d_hrmap;
  yoga_obj<int>            *d_binmap;
  yoga_obj<int>            *d_imamap;
  yoga_obj<float>          *d_offsets;
  yoga_obj<float>          *d_fluxPerSub;
  yoga_obj<float>          *d_pupil;
  yoga_obj<float>          *d_lgskern;
  yoga_obj<cuFloatComplex> *d_ftlgskern;
  yoga_obj<cuFloatComplex> *d_fttotim;

  yoga_obj<float>          *d_slopes;     
  yoga_obj<float>          *d_refslopes;

  // this is for centroiding
  yoga_obj<float>          *d_validpix;     
  yoga_obj<int>            *d_validindx;

  yoga_obj<float>          *d_weight; 
    
  yoga_obj<float>          *d_corrfct;     
  yoga_obj<cuFloatComplex> *d_corrfft1;     
  yoga_obj<cuFloatComplex> *d_corrfft2;     
  yoga_obj<float>          *d_corrnorm;     
  yoga_obj<float>          *d_corr;     

 public:
  yoga_wfs(long nxsub, long nvalid, long npix, long nphase, long nrebin, long nfft, long ntot, long npup,float pdiam,float nphotons, int lgs);
  yoga_wfs(const yoga_wfs& wfs);
  ~yoga_wfs();

  int wfs_initarrays(int *phasemap,int *hrmap,int *imamap,int *binmap,float *offsets, float *pupil, float *fluxPerSub);
  int wfs_initgs(float xpos,float ypos,float lambda, float mag, long size,float noise, long seed);
  int load_kernels(float *lgskern, int device);
  int load_corrfct(float *corrfct, int device);
  int init_nmax(int nmax);  
  int sensor_trace(yoga_atmos *yatmos);
  int comp_image(int device);
  int slopes_geom(int type);
  int get_cog();  
  int get_tcog(float threshold);  
  int get_bpcog(int nmax);
  int get_corr(int device);
  int get_nmax();
};

class yoga_sensors {
 public:
  int                 nsensors;
  vector<yoga_wfs *>  d_wfs;
     
 public:
  yoga_sensors(int nwfs,long *nxsub,long *nvalid,long *npix,long *nphase, long *nrebin,long *nfft, long *ntot, long npup, float *pdiam, float *nphot,int *lgs);
  ~yoga_sensors();

  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, float *noise, long *seed);
  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, float *noise);
  int sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size);
};

// General utilities
int fillcamplipup(cuFloatComplex *amplipup, float *phase,float *offset, float *mask, int *indx, int Nfft, int Npup, int Nsub, int npup, int device);
int fillbincube(float *bcube, float *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin, int Nsub, int device);
int fillbinimg(float *bimage, float *bcube, int *indxbin, int Npix, int Nsub, bool add,int device);
int indexfill(float *d_odata,float *d_idata,int *indx,int nx, int Nx,int N,int device);
int cfillrealp(cuFloatComplex *d_odata,float *d_idata,int N,int device);
int convolve(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device);
int cgetrealp(float *d_odata,cuFloatComplex *d_idata,int N,int device);
int subap_norm(float *d_odata,float *d_idata,float *fact,float *norm,float nphot,int n, int N,int device);
int fill_corr(cuFloatComplex *d_out, float *d_in, int npix_in, int npix_out, int N, int device);
int fillval_corr(cuFloatComplex *d_out, float val, int npix_in, int npix_out, int N, int device);
int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, int device);
int corr_norm(float *d_odata,float *d_idata,int N,int device);
int correl(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device);

// CUDA templates
template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata);
template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata, T thresh);
template <class T> void phase_reduce(int threads, int blocks, T *d_idata, T *d_odata, int *indx, T *offset, T alpha);
template <class T> void phase_derive(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, int *indx, T *offset, T *mask, T alpha, float *fluxPerSub);
template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, T *alpha);
template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, T *alpha, T thresh);
template <class T> void subap_sortmax(int size, int threads, int blocks, T *d_idata, T *d_odata, int *values, int nmax);
template <class T> void subap_centromax(int threads, int blocks, T *d_idata, T *d_odata, int npix, int nmax);
template <class T> int  roll(T *idata,int N,int M,int nim);

#endif // _YOGA_WFS_H_

