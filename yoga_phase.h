#ifndef _YOGA_PHASE_H_
#define _YOGA_PHASE_H_

// this is the generic class for a phase
// contains a yoga obj for the phase screen itself
// the screen size (assumed square)
// the following is only initialized on demand :
// an array of zernike coeffs on which the phase can be decomposed
// an array of zernike polynomials (yoga_object)
// a matrix to decompose the phase on zernike coeffs

#include <yoga.h>
#include <yoga_obj.h>

using namespace std;

class yoga_phase {
 public:

  yoga_obj<float>  *d_screen;     
  long             screen_size;
  float            *zerCoeff;
  yoga_obj<float>  *zernikes;     
  yoga_obj<float>  *mat;     

 public:
  yoga_phase(long size);
  yoga_phase(const yoga_phase& phase);
  ~yoga_phase();

};

template<class T> int launch_tcopy(T *d_odata,T *d_idata,int nx, int ny);

int phase_copy(yoga_phase *phase1,yoga_phase *phase2); 

#endif // _YOGA_PHASE_H_
