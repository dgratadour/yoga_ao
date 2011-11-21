#ifndef _YOGA_AO_UTILS_H_
#define _YOGA_AO_UTILS_H_
#include <yoga.h>

int cfillrealp(cuFloatComplex *d_odata,float *d_idata,int N,int device);
int cgetrealp(float *d_odata,cuFloatComplex *d_idata,int N,int device);
int abs2(float *d_odata, cuFloatComplex *d_idata, int N, int device);
int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, int device);
int subap_norm(float *d_odata,float *d_idata,float *fact,float *norm,float nphot,int n, int N,int device);
int fillindx(float *d_odata,float *d_idata,int *indx,int N, int device);
int fillarr2d(float *d_odata,float *d_idata,int x0, int Ncol,int NC, int N,int device);
int getarr2d(float *d_odata,float *d_idata,int x0, int Ncol,int NC, int N, int device);
int addai(float *d_odata,float *i_data, int i,int sgn, int N, int device);

// templates
template <class T> int  roll(T *idata,int N,int M,int nim);

#endif // _YOGA_AO_UTILS_H_

