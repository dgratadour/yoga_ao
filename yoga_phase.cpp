#include <yoga_phase.h>
#include "yapi.h"

yoga_phase::yoga_phase(long size)
{
  this->screen_size = size;
  
  long *dims_data2 = new long[3];
  dims_data2[0] = 2; dims_data2[1] = this->screen_size; 
  dims_data2[2] = this->screen_size; 
  
  this->d_screen = new yoga_obj<float>(dims_data2);
}

yoga_phase::~yoga_phase()
{
  delete this->d_screen;
}

int phase_copy(yoga_phase *phase1,yoga_phase *phase2) 
{

  launch_tcopy(phase1->d_screen->d_data,phase2->d_screen->d_data,
	       phase1->d_screen->dims_data[1],phase1->d_screen->dims_data[2]);

  return EXIT_SUCCESS;
}
