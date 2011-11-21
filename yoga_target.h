#ifndef _YOGA_TARGET_H_
#define _YOGA_TARGET_H_

#include <map>
#include <string>
#include <vector>
#include <yoga.h>
#include <yoga_turbu.h>
#include <yoga_lgs.h>

using namespace std;
using std::string;

typedef std::pair<std::string,float> type_screen;

class yoga_source {
 public:

  int                       device;      // device # 
  float                     tposx;       // x position of target on the sky  
  float                     tposy;       // y position of target on the sky
  long                      npos;        // number of points in the pupil
  float                     mag;         // brightness of target
  float                     lambda;      // imaging lambda
  float                     zp;          // imaging zero point
  bool                      lgs;         // flag for lgs
  string                    type;        // type of source : target / wfs
  yoga_phase                *d_phase;    // phase for this target
  yoga_lgs                  *d_lgs;      // the lgs object
  yoga_obj<float>           *object;     // the object intensity map
  yoga_obj<float>           *d_image;    // the resulting image for this target
  yoga_obj<cuFloatComplex>  *d_amplipup; // the complex amplitude in the pupil plane
  map<type_screen,float>    xoff;        // x reference for raytracing
  map<type_screen,float>    yoff;        // y reference for raytracing

 public:
  yoga_source(float xpos,float ypos,float lambda,float mag,long size,string type, int device);
  yoga_source(const yoga_source& source);
  ~yoga_source();
  int add_layer(string type,float alt,float xoff, float yoff);
  int remove_layer(string type,float alt);
  int get_phase(float *dest);
  int raytrace(yoga_atmos *atmos);
  int raytrace_shm(yoga_atmos *atmos);
  int comp_image(float *mask);
};

class yoga_target {
 public:
  int                    ntargets;
  vector<yoga_source *>  d_targets;
     
 public:
  yoga_target(int ntargets,float *xpos,float *ypos,float *lambda,float *mag,long *sizes, int device);
  ~yoga_target();

  int get_phase(int ntarget,float *dest);

};

int target_texraytrace(float *d_odata,float *d_idata,int nx, int ny,int Nx, int Ny, float xoff, 
		       float yoff, int Ntot, cudaChannelFormatDesc channelDesc, int device);
int target_raytrace(float *d_odata,float *d_idata,int nx, int ny,int Nx,float xoff, float yoff,  
		    int device);
int fillampli(cuFloatComplex *d_odata,float *d_idata, float *mask,int nx, int ny, int Nx, int device);
int fillpupil(cuFloatComplex *d_odata, float *mask,int nx, int ny, int Nx, int device);
int fft_goodsize(long size);

#endif // _YOGA_TARGET_H_
