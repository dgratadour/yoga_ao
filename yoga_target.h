#ifndef _YOGA_TARGET_H_
#define _YOGA_TARGET_H_

// this is the generic class for a phase
// contains a yoga obj for the phase screen itself
// the screen size (assumed square)
// the following is only initialized on demand :
// an array of zernike coeffs on which the phase can be decomposed
// an array of zernike polynomials (yoga_object)
// a matrix to decompose the phase on zernike coeffs

#include <map>
#include <string>
#include <vector>
#include <yoga.h>
#include <yoga_turbu.h>

using namespace std;
using std::string;

typedef std::pair<std::string,float> type_screen;

class yoga_source {
 public:

  float  tposx;                      // x position of target on the sky  
  float  tposy;                      // y position of target on the sky
  long   npos;                      // x position of target on the sky  
  map<type_screen,yoga_obj<float> *> rays_posx;  // light rays over the pupil
  map<type_screen,yoga_obj<float> *> rays_posy;  // light rays over the pupil
  yoga_phase *d_phase;               // phase for this target
  float  mag;                        // brightness of target
  float  lambda;                     // imaging lambda
  float  zp    ;                     // imaging zero point
  bool   lgs;                        // flag for lgs
  string type;                       // type of source : target / wfs
  yoga_obj<float> *object;           // the object intensity map
  yoga_obj<float> *d_image;          // the resulting image for this target
  yoga_obj<cuFloatComplex> *d_amplipup;       // the complex amplitude in the pupil plane

 public:
  yoga_source(float xpos,float ypos,float lambda,float mag,long size,string type);
  yoga_source(const yoga_source& source);
  ~yoga_source();
  int add_layer(string type,float alt,float *xref,float *yref,long size);
  int remove_layer(string type,float alt);
  int get_phase(float *dest);
  int raytrace(yoga_atmos *atmos);
  int comp_image(float *mask);
};

class yoga_target {
 public:
  int                    ntargets;
  vector<yoga_source *>  d_targets;
     
 public:
  yoga_target(int ntargets,float *xpos,float *ypos,float *lambda,float *mag,long *sizes);
  ~yoga_target();

  int set_layer(int ntarget,string type,float alt,float *xref,float *yref,long size);
  int get_phase(int ntarget,float *dest);

};

int launch_raytrace(float *d_odata,float *d_idata,int nx, int ny,float *xref, float *yref, int Nx);
int launch_fillampli(cuFloatComplex *d_odata,float *d_idata, float *mask,int nx, int ny, int Nx);
int launch_fillpupil(cuFloatComplex *d_odata, float *mask,int nx, int ny, int Nx);
int launch_abs2(float *d_odata, cuFloatComplex *d_idata, int nx, int ny);
int fft_goodsize(long size);

#endif // _YOGA_TARGET_H_
