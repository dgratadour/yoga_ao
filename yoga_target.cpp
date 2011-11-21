#include <yoga_ao_utils.h>
#include <yoga_target.h>
#include <yoga_phase.h>

int fft_goodsize(long size)
{
  int mradix = 2;
  float tmpf = 0.;
  long tmpl = 0;

  tmpf = logf(size)/logf(3);
  tmpl = (long)tmpf;
  tmpf -= tmpl;
  mradix = (tmpf > (logf(size)/logf(mradix) - (long)(logf(size)/logf(mradix))) ? 3 : mradix);

  tmpf = logf(size)/logf(5);
  tmpl = (long)tmpf;
  tmpf -= tmpl;
  mradix = (tmpf > (logf(size)/logf(mradix) - (long)(logf(size)/logf(mradix))) ? 5 : mradix);

  tmpf = logf(size)/logf(7);
  tmpl = (long)tmpf;
  tmpf -= tmpl;
  mradix = (tmpf > (logf(size)/logf(mradix) - (long)(logf(size)/logf(mradix))) ? 7 : mradix);

  return mradix;
}


yoga_source::yoga_source(float xpos,float ypos,float lambda, float mag, long size, string type, int device)
{
  tposx         = xpos;    
  tposy         = ypos;   
  this->lambda  = lambda;
  this->mag     = mag;
  this->npos    = size;
  this->d_phase = new yoga_phase(size);
  this->type    = type;
  this->device  = device;

  if (type != "wfs") {
    int mradix = fft_goodsize(size);
    
    int fft_size = pow(mradix,(long)(logf(size)/logf(mradix))+1);
    
    long *dims_data2 = new long[3];
    dims_data2[0] = 2; 
    dims_data2[1] = fft_size;
    dims_data2[2] = fft_size; 
    
    this->d_image = new yoga_obj<float>(dims_data2);
    this->d_amplipup = new yoga_obj<cuFloatComplex>(dims_data2);

    cufftSafeCall(cufftPlan2d(&(this->d_amplipup->plan), this->d_amplipup->dims_data[1],
			       this->d_amplipup->dims_data[2],CUFFT_C2C));
    this->d_amplipup->fft_on = true;
  }
}

yoga_source::~yoga_source()
{
  delete this->d_phase;
  if (this->type != "wfs") {
    delete this->d_image;
    delete this->d_amplipup;
  }
  this->xoff.clear();
  this->yoff.clear();
}

int yoga_source::add_layer(string type,float alt,float mxoff, float myoff)
{
    xoff[make_pair(type,alt)] = mxoff;
    yoff[make_pair(type,alt)] = myoff;

  return EXIT_SUCCESS;
}

int yoga_source::remove_layer(string type,float alt)
{
  xoff.erase(make_pair(type,alt));
  yoff.erase(make_pair(type,alt));
  
  return EXIT_SUCCESS;
}

int yoga_source::raytrace_shm(yoga_atmos *yatmos)
{
  cutilSafeCall(cudaMemset(this->d_phase->d_screen->d_data, 0, 
			   sizeof(float)*this->d_phase->d_screen->nb_elem));

  map<type_screen,float>::iterator p;
  p = xoff.begin();
  while (p != xoff.end()) {
    string types = p->first.first;
    if(types.find("atmos")==0) {
      float alt = p->first.second;
      map<float,yoga_tscreen *>::iterator ps;
      ps = yatmos->d_screens.find(alt);
      if (ps != yatmos->d_screens.end()) {
	target_texraytrace(this->d_phase->d_screen->d_data,
			   ps->second->d_tscreen->d_screen->d_data,
			   (int)d_phase->d_screen->dims_data[1],
			   (int)d_phase->d_screen->dims_data[2],
			   (int)ps->second->d_tscreen->d_screen->dims_data[1],
			   (int)ps->second->d_tscreen->d_screen->dims_data[2],
			   xoff[make_pair("atmos",alt)],
			   yoff[make_pair("atmos",alt)],
			   ps->second->d_tscreen->d_screen->nb_elem,
			   ps->second->channelDesc,
			   ps->second->device);
      }
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int yoga_source::raytrace(yoga_atmos *yatmos)
{
  cutilSafeCall(cudaMemset(this->d_phase->d_screen->d_data, 0, 
			   sizeof(float)*this->d_phase->d_screen->nb_elem));

  map<type_screen,float>::iterator p;
  p = xoff.begin();
  while (p != xoff.end()) {
    string types = p->first.first;
    if(types.find("atmos")==0) {
      float alt = p->first.second;
      map<float,yoga_tscreen *>::iterator ps;
      ps = yatmos->d_screens.find(alt);
      if (ps != yatmos->d_screens.end()) {
	target_raytrace(this->d_phase->d_screen->d_data,
			ps->second->d_tscreen->d_screen->d_data,
			(int)d_phase->d_screen->dims_data[1],
			(int)d_phase->d_screen->dims_data[2],
			(int)ps->second->d_tscreen->d_screen->dims_data[1],
			xoff[make_pair("atmos",alt)],
			yoff[make_pair("atmos",alt)],
			ps->second->device);
      }
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int yoga_source::comp_image(float *mask)
{
  cutilSafeCall(cudaMemset(this->d_amplipup->d_data, 0,sizeof(cuFloatComplex)*this->d_amplipup->nb_elem));
  
  fillampli(this->d_amplipup->d_data,this->d_phase->d_screen->d_data, mask,this->d_phase->d_screen->dims_data[1], 
	    this->d_phase->d_screen->dims_data[2],this->d_amplipup->dims_data[1],this->device);
  
  /*
  fillpupil(this->d_amplipup->d_data, mask, this->d_phase->d_screen->dims_data[1], 
		   this->d_phase->d_screen->dims_data[2], this->d_amplipup->dims_data[1],this->device);
  */
  
  yoga_fft(this->d_amplipup->d_data,this->d_amplipup->d_data,1,this->d_amplipup->plan);
  
  
  abs2(this->d_image->d_data,this->d_amplipup->d_data,this->d_image->dims_data[1]*this->d_image->dims_data[2],this->device);
  
  return EXIT_SUCCESS;
}


yoga_target::yoga_target(int ntargets,float *xpos,float *ypos,float *lambda,float *mag, long *sizes, int device)
{
  this->ntargets = ntargets;

  for (int i=0;i<ntargets;i++) {
    d_targets.push_back(new yoga_source(xpos[i],ypos[i],lambda[i],mag[i],sizes[i],"target",device));
  }
}


yoga_target::~yoga_target()
{
  for (size_t idx = 0; idx < (this->d_targets).size(); idx++) {
   delete (this->d_targets)[idx];
  } 
}

