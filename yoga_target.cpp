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


yoga_source::yoga_source(float xpos,float ypos,float lambda, float mag, long size, string type)
{
  tposx         = xpos;    
  tposy         = ypos;   
  this->lambda  = lambda;
  this->mag     = mag;
  this->npos    = size;
  this->d_phase = new yoga_phase(size);
  this->type    = type;

  if (type != "wfs") {
    int mradix = fft_goodsize(size);
    
    int fft_size = pow(mradix,(long)(logf(size)/logf(mradix))+1);
    
    long *dims_data2 = new long[3];
    dims_data2[0] = 2; dims_data2[1] = fft_size;dims_data2[2] = fft_size; 
    
    this->d_image = new yoga_obj<float>(dims_data2);
    this->d_amplipup = new yoga_obj<cuFloatComplex>(dims_data2);
    cufftSafeCall( cufftPlan2d(&(this->d_amplipup->plan), this->d_amplipup->dims_data[1],
			       this->d_amplipup->dims_data[2],CUFFT_C2C));
    this->d_amplipup->fft_on = true;
  }
}

yoga_source::~yoga_source()
{
  delete this->d_phase;
  for ( map<type_screen,yoga_obj<float> *>::iterator it = rays_posx.begin(); it != rays_posx.end(); ++it )
    {
      delete it->second;
      it->second = 0;
    }  
  for ( map<type_screen,yoga_obj<float> *>::iterator it = rays_posy.begin(); it != rays_posy.end(); ++it )
    {
      delete it->second;
      it->second = 0;
    }  
  if (this->type != "wfs") {
    delete this->d_image;
    delete this->d_amplipup;
  }
}

int yoga_source::add_layer(string type,float alt,float *xref,float *yref,long size)
{
  if (size == this->npos) {
    long dims_data[2];
    dims_data[0]=1;
    dims_data[1] = size;
    rays_posx[make_pair(type,alt)] = new yoga_obj<float>(dims_data);
    (rays_posx[make_pair(type,alt)])->host2device(xref);
    rays_posy[make_pair(type,alt)] = new yoga_obj<float>(dims_data);
    rays_posy[make_pair(type,alt)]->host2device(yref);
  } else cout << "Wrong rays dimension" << endl;

  return EXIT_SUCCESS;
}

int yoga_source::remove_layer(string type,float alt)
{
 rays_posx.erase(make_pair(type,alt));
 rays_posy.erase(make_pair(type,alt));
  
  return EXIT_SUCCESS;
}

int yoga_source::raytrace(yoga_atmos *yatmos)
{
  cutilSafeCall(cudaMemset(this->d_phase->d_screen->d_data, 0, 
			   sizeof(float)*this->d_phase->d_screen->nb_elem));

  map<type_screen,yoga_obj<float> *>::iterator p;
  p = rays_posx.begin();
  while (p != rays_posx.end()) {
    string types = p->first.first;
    if(types.find("atmos")==0) {
      float alt = p->first.second;
      map<float,yoga_tscreen *>::iterator ps;
      ps = yatmos->d_screens.find(alt);
   
      if (ps != yatmos->d_screens.end())
	target_raytrace(this->d_phase->d_screen->d_data,
			ps->second->d_tscreen->d_screen->d_data,
			(int)d_phase->d_screen->dims_data[1],
			(int)d_phase->d_screen->dims_data[2],
			p->second->d_data,
			rays_posy[make_pair("atmos",alt)]->d_data,
			(int)ps->second->d_tscreen->d_screen->dims_data[1]);
    }
    p++;
  }

  return EXIT_SUCCESS;
}

int yoga_source::comp_image(float *mask)
{
  cutilSafeCall(cudaMemset(this->d_amplipup->d_data, 0,sizeof(cuFloatComplex)*this->d_amplipup->nb_elem));
  
  fillampli(this->d_amplipup->d_data,this->d_phase->d_screen->d_data, mask,this->d_phase->d_screen->dims_data[1], 
		   this->d_phase->d_screen->dims_data[2],this->d_amplipup->dims_data[1]);
  
  /*
  fillpupil(this->d_amplipup->d_data, mask, this->d_phase->d_screen->dims_data[1], 
		   this->d_phase->d_screen->dims_data[2], this->d_amplipup->dims_data[1]);
  */
  
  yoga_fft(this->d_amplipup->d_data,this->d_amplipup->d_data,1,this->d_amplipup->plan);
  
  
  abs2(this->d_image->d_data,this->d_amplipup->d_data,this->d_image->dims_data[1],
	      this->d_image->dims_data[2]);
  
  return EXIT_SUCCESS;
}


yoga_target::yoga_target(int ntargets,float *xpos,float *ypos,float *lambda,float *mag, long *sizes)
{
  this->ntargets = ntargets;

  for (int i=0;i<ntargets;i++) {
    d_targets.push_back(new yoga_source(xpos[i],ypos[i],lambda[i],mag[i],sizes[i],"target"));
  }
}


yoga_target::~yoga_target()
{
  for (size_t idx = 0; idx < (this->d_targets).size(); idx++) {
   delete (this->d_targets)[idx];
  } 
}

