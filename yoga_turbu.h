#ifndef _YOGA_TURBU_H_
#define _YOGA_TURBU_H_

#include <vector>
#include <map>
#include <yoga_phase.h>

using namespace std;

class yoga_tscreen {
 public:

  yoga_phase             *d_tscreen;     
  yoga_obj<float>        *d_tscreen_o;     
  yoga_obj<float>        *d_A;
  yoga_obj<float>        *d_B;
  yoga_obj<unsigned int> *d_istencilx;
  yoga_obj<unsigned int> *d_istencily;
  yoga_obj<float>        *d_z;
  yoga_obj<float>        *d_noise;
  yoga_obj<float>        *d_ytmp;     
  long                   screen_size;
  float                  amplitude;
  float                  altitude;
  float                  windspeed;
  float                  winddir;
  float                  deltax;
  float                  deltay;
  float                  accumx;
  float                  accumy;

 public:
  yoga_tscreen(long size, long size2, float amplitude, float altitude, float windspeed, float winddir, float deltax, float deltay);
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
  yoga_atmos(int nscreens,float *r0,long *size, long *size2, float *altitude, float *windspeed, float *winddir, float *deltax, float *deltay, float *pupil);
  ~yoga_atmos();

  int init_screen(float alt,float *h_A, float *h_B, unsigned int *h_istencilx,unsigned int *h_istencily, int seed);
  int move();
};

extern "C" {
};

#endif // _YOGA_TURBU_H_

