#ifndef _YOGA_TURBU_H_
#define _YOGA_TURBU_H_

#include <vector>
#include <map>
#include <yoga_ao_utils.h>
#include <yoga_phase.h>

using namespace std;

class yoga_tscreen {
 public:

  int                    device;       // The device # 
  yoga_phase             *d_tscreen;   // The phase screen   
  yoga_obj<float>        *d_tscreen_o; // Additional space of the same size as the phase screen
  yoga_obj<float>        *d_A;         // A matrix for extrusion
  yoga_obj<float>        *d_B;         // B matrix for extrusion
  yoga_obj<unsigned int> *d_istencilx; // stencil for column extrusion
  yoga_obj<unsigned int> *d_istencily; // stencil for row extrusion
  yoga_obj<float>        *d_z;         // tmp array for extrusion process
  yoga_obj<float>        *d_noise;     // tmp array containing random numbers
  yoga_obj<float>        *d_ytmp;      // contains the extrude update (row or column)
  long                   screen_size;  // size of phase screens
  float                  amplitude;    // amplitude for extrusion (r0^-5/6)
  float                  altitude;     
  float                  windspeed;
  float                  winddir;
  float                  deltax;       // number of row to extrude per iteration
  float                  deltay;       // number of row to extrude per iteration
  //internal
  float                  accumx;
  float                  accumy;
  cudaChannelFormatDesc  channelDesc;  // Channel descriptor for texture memory access

 public:
  yoga_tscreen(long size, long size2, float amplitude, float altitude, float windspeed, float winddir, float deltax, float deltay,int device);
  yoga_tscreen(const yoga_tscreen& tscreen);
  ~yoga_tscreen();

  int init_screen(float *h_A, float *h_B, unsigned int *h_istencilx,unsigned int *h_istencily, int seed);
  int extrude(int dir);
};


class yoga_atmos {
 public:
  int                    nscreens;
  map<float,yoga_tscreen *> d_screens;
  float                  r0;
  yoga_obj<float>        *d_pupil;     
     
 public:
  yoga_atmos(int nscreens,float *r0,long *size, long *size2, float *altitude, float *windspeed, float *winddir, float *deltax, float *deltay, float *pupil, int device);
  ~yoga_atmos();

  int init_screen(float alt,float *h_A, float *h_B, unsigned int *h_istencilx,unsigned int *h_istencily, int seed);
  int move();
};

extern "C" {
};

#endif // _YOGA_TURBU_H_

