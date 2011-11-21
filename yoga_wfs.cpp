#include <yoga_wfs.h>
#include <yoga_ao_utils.h>

/*
  Algorithm for sh wfs :
  copy (phase+offset)*mask into camplipup
  do multi-fft on camplipup to ctubrim
  if larger fov : fill larger image
  fill bimage with noise
  add binned pixels from cturbim to bimage
 */
yoga_wfs::yoga_wfs(long nxsub, long nvalid, long npix, long nphase, long nrebin, long nfft, 
		   long ntot, long npup,float pdiam,float nphotons, int lgs, int device)
{
  this->nxsub   = nxsub;
  this->nvalid  = nvalid;
  this->npix    = npix;
  this->nphase  = nphase;
  this->nrebin  = nrebin;
  this->nfft    = nfft;
  this->ntot    = ntot;
  this->npup    = npup;
  this->subapd  = pdiam;
  this->nphot   = nphotons;
  this->lgs     = (lgs == 1 ? true : false);
  this->device  = device;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2; 
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  if (ntot != nfft) {
    dims_data3[1] = ntot; dims_data3[2] = ntot; dims_data3[3] = nvalid;  
    this->d_totimg = new yoga_obj<float>(dims_data3);
  } 

  dims_data2[1] = npix*nxsub; dims_data2[2] =npix*nxsub ; 
  this->d_binimg = new yoga_obj<float>(dims_data2);

  dims_data3[1] = nfft; dims_data3[2] =nfft; dims_data3[3] = nvalid;  
  this->d_hrimg = new yoga_obj<float>(dims_data3);
  this->d_camplipup = new yoga_obj<cuFloatComplex>(dims_data3);
  this->d_camplifoc = new yoga_obj<cuFloatComplex>(dims_data3);

  int mdims[2];

  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];

  cufftSafeCall( cufftPlanMany(&(this->d_camplipup->plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
			       (int)dims_data3[3]));
  this->d_camplipup->fft_on = true;

  dims_data3[1] = npix; dims_data3[2] = npix; dims_data3[3] = nvalid;  
  this->d_bincube = new yoga_obj<float>(dims_data3);

  if (this->lgs) {
    dims_data3[1] = ntot; dims_data3[2] =ntot; dims_data3[3] = nvalid;  
    this->d_fttotim   = new yoga_obj<cuFloatComplex>(dims_data3);
    mdims[0] = (int)dims_data3[1];
    mdims[1] = (int)dims_data3[2];
    cufftSafeCall( cufftPlanMany(&(this->d_fttotim->plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
				 (int)dims_data3[3]));
    this->d_fttotim->fft_on = true;
  }

  dims_data2[1] = nphase*nphase; dims_data2[2] = nvalid; 
  this->d_phasemap = new yoga_obj<int>(dims_data2);

  if (this->ntot != this->nfft) {
    dims_data1[1] = nfft * nfft;  
    this->d_hrmap = new yoga_obj<int>(dims_data1);
  }

  dims_data2[1] = npup; dims_data2[2] = npup; 
  this->d_pupil = new yoga_obj<float>(dims_data2);

  dims_data2[1] = nphase; dims_data2[2] = nphase; 
  this->d_offsets = new yoga_obj<float>(dims_data2);

  dims_data2[1] = nxsub; dims_data2[2] = nxsub; 
  this->d_isvalid = new yoga_obj<int>(dims_data2);

  dims_data1[1] = nvalid; 
  this->d_validsubsx = new yoga_obj<int>(dims_data1);
  this->d_validsubsy = new yoga_obj<int>(dims_data1);

  dims_data1[1] = nxsub;
  this->d_istart = new yoga_obj<int>(dims_data1);
  this->d_jstart = new yoga_obj<int>(dims_data1);

  dims_data2[1] = nrebin*nrebin; dims_data2[2] = npix * npix; 
  this->d_binmap = new yoga_obj<int>(dims_data2);

  dims_data1[1] = nvalid;
  this->d_subsum = new yoga_obj<float>(dims_data1);
  this->d_fluxPerSub = new yoga_obj<float>(dims_data1);

  dims_data1[1] = 2*nvalid;
  this->d_slopes = new yoga_obj<float>(dims_data1);

  this->d_validpix = 0L;
  this->d_validindx = 0L;

  this->d_weight = 0L;

  this->d_corrfct  = 0L;
  this->d_corrfft1 = 0L;
  this->d_corrfft2 = 0L;
  this->d_corrnorm = 0L;
  this->d_corr     = 0L;


}

yoga_wfs::~yoga_wfs()
{
  delete this->d_binimg;
  if (this->ntot != this->nfft) {
    delete this->d_totimg;
    delete this->d_hrmap;
  }
  delete this->d_hrimg;
  delete this->d_camplipup;
  delete this->d_camplifoc;
  delete this->d_binmap;
  delete this->d_phasemap;
  delete this->d_offsets;
  delete this->d_pupil;
  delete this->d_subsum;
  delete this->d_fluxPerSub;
  delete this->d_slopes;
  delete this->d_gs;
  delete this->d_isvalid;
  delete this->d_validsubsx;
  delete this->d_validsubsy;
  delete this->d_istart;
  delete this->d_jstart;

  if (this->d_validpix != 0L) delete this->d_validpix;
  if (this->d_validindx!= 0L) delete this->d_validindx;

  if (this->d_weight != 0L) delete this->d_weight;

  if (this->d_corrfct != 0L) delete this->d_corrfct;
  if (this->d_corrfft1 != 0L) delete this->d_corrfft1;
  if (this->d_corrfft2 != 0L) delete this->d_corrfft2;
  if (this->d_corrnorm != 0L) delete this->d_corrnorm;
  if (this->d_corr != 0L) delete this->d_corr;

}

int yoga_wfs::wfs_initgs(float xpos,float ypos,float lambda, float mag, long size,float noise,long seed)
{
  this->d_gs = new yoga_source(xpos,ypos,lambda,mag,size,"wfs",this->device);
  this->noise = noise;
  if (noise > 0) {
    this->d_binimg->init_prng(seed);
    this->d_binimg->prng('N',noise);
  }
  if (this->lgs) {
    this->d_gs->d_lgs  = new yoga_lgs(this->nvalid,this->ntot);
    this->d_gs->lgs = this->lgs;
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::wfs_initarrays(int *phasemap,int *hrmap, int *binmap,float *offsets, 
			     float *pupil, float *fluxPerSub, int *isvalid, int *validsubsx, int *validsubsy, 
			     int *istart, int *jstart)
{
  this->d_phasemap->host2device(phasemap);
  this->d_offsets->host2device(offsets);
  this->d_pupil->host2device(pupil);
  this->d_binmap->host2device(binmap);
  this->d_fluxPerSub->host2device(fluxPerSub);
  if (this->ntot != this->nfft) this->d_hrmap->host2device(hrmap);
  this->d_validsubsx->host2device(validsubsx);
  this->d_validsubsy->host2device(validsubsy);
  this->d_isvalid->host2device(isvalid);
  this->d_istart->host2device(istart);
  this->d_jstart->host2device(jstart);
  return EXIT_SUCCESS;
}

int yoga_wfs::load_kernels(float *lgskern)
{
  if (this->lgs) this->d_gs->d_lgs->load_kernels(lgskern,this->device);

  return EXIT_SUCCESS;
}

int yoga_wfs::load_corrfct(float *corrfct)
{
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;
  dims_data3[1] = npix; dims_data3[2] = npix; dims_data3[3] = nvalid;  
  if (this->d_corrfct  == 0L) this->d_corrfct = new yoga_obj<float>(dims_data3);

  dims_data3[1] = 2*npix; dims_data3[2] = 2*npix; 
  int mdims[2];
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];

  if (this->d_corrfft1 == 0L) {
    this->d_corrfft1 = new yoga_obj<cuFloatComplex>(dims_data3);
    cufftSafeCall(cufftPlanMany(&(this->d_corrfft1->plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
				(int)dims_data3[3]));
    this->d_corrfft1->fft_on = true;
  }

  if (this->d_corrfft2 == 0L) this->d_corrfft2 = new yoga_obj<cuFloatComplex>(dims_data3);

  dims_data3[1] = 2*npix-1; dims_data3[2] = 2*npix-1;
  if (this->d_corrnorm  == 0L) this->d_corrnorm = new yoga_obj<float>(dims_data3);
  if (this->d_corr  == 0L) this->d_corr = new yoga_obj<float>(dims_data3);

  this->d_corrfct->host2device(corrfct);

  fill_corr(this->d_corrfft2->d_data,this->d_corrfct->d_data,this->npix*this->npix,
	    this->npix*this->npix*4,this->d_corrfct->nb_elem,this->device);

  yoga_fft(this->d_corrfft2->d_data,this->d_corrfft2->d_data,1,this->d_corrfft1->plan);

  // here we compute the normalization map
  fillval_corr(this->d_corrfft1->d_data,1.0f,this->npix*this->npix,this->npix*this->npix*4,
	       this->d_corrfct->nb_elem,this->device);

  yoga_fft(this->d_corrfft1->d_data,this->d_corrfft1->d_data,1,this->d_corrfft1->plan);

  abs2c(this->d_corrfft1->d_data,this->d_corrfft1->d_data ,this->d_corrfft1->nb_elem ,this->device);

  yoga_fft(this->d_corrfft1->d_data,this->d_corrfft1->d_data,-1,this->d_corrfft1->plan);

  roll(this->d_corrfft1->d_data,this->d_corrfft1->dims_data[1],this->d_corrfft1->dims_data[2],
       this->d_corrfft1->dims_data[3]);

  cgetrealp(this->d_corrnorm->d_data,&(this->d_corrfft1->d_data[this->d_corrfft1->dims_data[1]+1]),
	    this->d_corrnorm->nb_elem,this->device);

  return EXIT_SUCCESS;
}

int yoga_wfs::init_nmax(int nmax)
{
  if (this->d_validpix != 0L) delete this->d_validpix;
  if (this->d_validindx!= 0L)  delete  this->d_validindx;

  long *dims_data2 = new long[3];
  dims_data2[0] = 2; 
  dims_data2[1] = nmax; dims_data2[2] = this->nvalid; 

  this->d_validpix = new yoga_obj<float>(dims_data2);
  this->d_validindx = new yoga_obj<int>(dims_data2);

  return EXIT_SUCCESS;
}

int yoga_wfs::sensor_trace(yoga_atmos *yatmos)
{
  //do raytracing to get the phase
  this->d_gs->raytrace(yatmos);

  return EXIT_SUCCESS;
}

int yoga_wfs::comp_image()
{
  //fill cube of complex ampli with exp(i*phase)
  /*
  fillcamplipup(this->d_camplifoc->d_data,this->d_gs->d_phase->d_screen->d_data, 
		this->d_offsets->d_data,this->d_pupil->d_data, this->d_phasemap->d_data, 
		this->nfft,this->nphase * this->nphase,this->nvalid,this->nphase,this->device);
  */
  fillcamplipup2(this->d_camplipup->d_data,this->d_gs->d_phase->d_screen->d_data, 
		 this->d_offsets->d_data,this->d_pupil->d_data, this->d_istart->d_data, this->d_jstart->d_data, 
		 this->d_validsubsx->d_data,this->d_validsubsy->d_data,this->nphase,
		 this->d_gs->d_phase->d_screen->dims_data[1],this->nfft,this->nphase*this->nphase*this->nvalid,
		 this->device);

  //do fft of the cube  
  yoga_fft(this->d_camplipup->d_data,this->d_camplifoc->d_data,1,this->d_camplipup->plan);

  //get the hrimage by taking the | |^2
  abs2(this->d_hrimg->d_data,this->d_camplifoc->d_data,this->d_hrimg->dims_data[1] 
       * this->d_hrimg->dims_data[2] * this->d_hrimg->dims_data[3], this->device);

  //set bincube to 0
  cutilSafeCall(cudaMemset(this->d_bincube->d_data, 0,sizeof(float)*this->d_bincube->nb_elem));
  
  // increase fov if required
  // and fill bincube with data from hrimg
  if (this->ntot != this->nfft) {
    
    cutilSafeCall(cudaMemset(this->d_totimg->d_data, 0, 
			     sizeof(float)*this->d_totimg->nb_elem));
    
    indexfill(this->d_totimg->d_data,this->d_hrimg->d_data,this->d_hrmap->d_data,
	      this->nfft,this->ntot,this->d_hrimg->nb_elem,this->device);

    if (this->lgs) {
      cutilSafeCall(cudaMemset(this->d_fttotim->d_data, 0, 
			       sizeof(cuFloatComplex)*this->d_fttotim->nb_elem));

      cfillrealp(this->d_fttotim->d_data,this->d_totimg->d_data,this->d_totimg->nb_elem,this->device);

      yoga_fft(this->d_fttotim->d_data,this->d_fttotim->d_data,1,this->d_fttotim->plan);

      convolve(this->d_fttotim->d_data,this->d_gs->d_lgs->d_ftlgskern->d_data,this->d_fttotim->nb_elem,this->device);

      yoga_fft(this->d_fttotim->d_data,this->d_fttotim->d_data,-1,this->d_fttotim->plan);

      cgetrealp(this->d_totimg->d_data,this->d_fttotim->d_data,this->d_fttotim->nb_elem,this->device);
      // note : i'm loosing time here need to rewrite fillbincube ...
    }

    fillbincube(this->d_bincube->d_data,this->d_totimg->d_data,this->d_binmap->d_data, 
		this->ntot * this->ntot,this->npix * this->npix, this->nrebin * this->nrebin,
		this->nvalid,this->device);
  } else {
    if (this->lgs) {
      cutilSafeCall(cudaMemset(this->d_fttotim->d_data, 0, 
			       sizeof(cuFloatComplex)*this->d_fttotim->nb_elem));

      cfillrealp(this->d_fttotim->d_data,this->d_hrimg->d_data,this->d_hrimg->nb_elem,this->device);

      yoga_fft(this->d_fttotim->d_data,this->d_fttotim->d_data,1,this->d_fttotim->plan);

      convolve(this->d_fttotim->d_data,this->d_gs->d_lgs->d_ftlgskern->d_data,this->d_fttotim->nb_elem,this->device);

      yoga_fft(this->d_fttotim->d_data,this->d_fttotim->d_data,-1,this->d_fttotim->plan);

      cgetrealp(this->d_hrimg->d_data,this->d_fttotim->d_data,this->d_fttotim->nb_elem,this->device);
      // note : i'm loosing time here need to rewrite fillbincube ...
    }

    fillbincube(this->d_bincube->d_data,this->d_hrimg->d_data ,this->d_binmap->d_data, 
		this->nfft * this->nfft,this->npix * this->npix, this->nrebin * this->nrebin,
		this->nvalid,this->device);
  }

  // normalize images :
  // get the sum value per subap
  subap_reduce(this->d_bincube->nb_elem,this->npix * this->npix,this->nvalid,this->d_bincube->d_data,
	       this->d_subsum->d_data);
  // multiply each subap by nphot*fluxPersub/sumPerSub
  subap_norm(this->d_bincube->d_data,this->d_bincube->d_data,this->d_fluxPerSub->d_data,
	     this->d_subsum->d_data,this->nphot,this->npix * this->npix,this->d_bincube->nb_elem ,this->device);

  //set binned image to noise   
  if (this->noise > 0) {
    this->d_binimg->prng('N',this->noise);
    fillbinimg(this->d_binimg->d_data,this->d_bincube->d_data,this->npix,this->nvalid,this->npix*this->nxsub,
	       this->d_validsubsx->d_data,this->d_validsubsy->d_data,true,this->device);
  } else {
    //fill binned image with data from bincube
    fillbinimg(this->d_binimg->d_data,this->d_bincube->d_data,this->npix,this->nvalid,this->npix*this->nxsub,
	       this->d_validsubsx->d_data,this->d_validsubsy->d_data,false,this->device);
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::slopes_geom(int type)
{
/*
    // computing subaperture phase difference at edges

 todo : integrale( x * phase ) / integrale (x^2);

*/
  if (type == 0) {
    float alpha = this->d_gs->lambda / this->subapd  / 4.84814 / 2 / 3.14159265
      * (float)this->nphase / (float)(this->nphase - 1);
    phase_reduce(this->nphase,this->nvalid,this->d_gs->d_phase->d_screen->d_data,
		 this->d_slopes->d_data,this->d_phasemap->d_data,alpha);
  }

  if (type == 1) {
    float alpha = this->d_gs->lambda / this->subapd  / 4.84814 / 2 / 3.14159265
      * (float)this->nphase / (float)(this->nphase - 1);
    phase_derive(this->nphase * this->nphase * this->nvalid,this->nphase * this->nphase,
		 this->nvalid,this->nphase,this->d_gs->d_phase->d_screen->d_data,
		 this->d_slopes->d_data,this->d_phasemap->d_data,this->d_pupil->d_data,
		 alpha,this->d_fluxPerSub->d_data);
  }

  return EXIT_SUCCESS;
}

int yoga_wfs::get_cog()
{
  // simple cog
  subap_reduce(this->d_bincube->nb_elem,this->npix * this->npix,this->nvalid,this->d_bincube->d_data,
	       this->d_subsum->d_data);
  
  get_centroids(this->d_bincube->nb_elem,this->npix * this->npix,this->nvalid,this->npix,
		this->d_bincube->d_data,this->d_slopes->d_data,this->d_subsum->d_data);

  return EXIT_SUCCESS;
}

int yoga_wfs::get_tcog(float threshold)
{
  // thresholded cog
  // in this case pixels < threshold are put to 0
  subap_reduce(this->d_bincube->nb_elem,this->npix * this->npix,this->nvalid,this->d_bincube->d_data,
	       this->d_subsum->d_data,threshold);
  
  get_centroids(this->d_bincube->nb_elem,this->npix * this->npix,this->nvalid,this->npix,
		this->d_bincube->d_data,this->d_slopes->d_data,this->d_subsum->d_data,threshold);

  return EXIT_SUCCESS;
}

int yoga_wfs::get_bpcog(int nmax)
{
  // brightest pixels cog
  subap_centromax(this->npix * this->npix,this->nvalid,this->d_bincube->d_data,this->d_slopes->d_data,
		  this->npix,nmax);
  return EXIT_SUCCESS;
}

int yoga_wfs::get_nmax()
{

  subap_sortmax(this->d_bincube->nb_elem,this->npix * this->npix ,this->nvalid,this->d_bincube->d_data, 
		this->d_validpix->d_data,this->d_validindx->d_data ,this->d_validpix->dims_data[1] );

  return EXIT_SUCCESS;
}

int yoga_wfs::get_corr()
{
  // correlation algorithm
  fill_corr(this->d_corrfft1->d_data,this->d_bincube->d_data,this->npix*this->npix,
	    this->npix*this->npix*4,this->d_bincube->nb_elem,this->device);

  yoga_fft(this->d_corrfft1->d_data,this->d_corrfft1->d_data,1,this->d_corrfft1->plan);

  correl(this->d_corrfft1->d_data,this->d_corrfft2->d_data,this->d_corrfft1->nb_elem, this->device);

  yoga_fft(this->d_corrfft1->d_data,this->d_corrfft1->d_data,-1,this->d_corrfft1->plan);

  roll(this->d_corrfft1->d_data,this->d_corrfft1->dims_data[1],this->d_corrfft1->dims_data[2],
       this->d_corrfft1->dims_data[3]);

  cgetrealp(this->d_corr->d_data,&(this->d_corrfft1->d_data[this->d_corrfft1->dims_data[1]+1]),
	    this->d_corr->nb_elem,this->device);

  //here need to normalize
  corr_norm(this->d_corr->d_data,this->d_corrnorm->d_data,this->d_corr->nb_elem,this->device);

  // need to find max for each subap

  // do parabolic interpolation


  return EXIT_SUCCESS;
}

yoga_sensors::yoga_sensors(int nwfs,long *nxsub,long *nvalid,long *npix,long *nphase, long *nrebin,
			   long *nfft, long *ntot,long npup,float *pdiam, float *nphot,  int *lgs, int device)
{
  this->nsensors = nwfs;

  for (int i=0;i<nwfs;i++) {
    d_wfs.push_back(new yoga_wfs(nxsub[i],nvalid[i],npix[i],nphase[i],nrebin[i],nfft[i],ntot[i],npup,
				 pdiam[i],nphot[i],lgs[i],device));
  }
}

yoga_sensors::~yoga_sensors()
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
   delete (this->d_wfs)[idx];
  } 
}

int yoga_sensors::sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, 
				 float *noise, long *seed)
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx],ypos[idx],lambda[idx],mag[idx],size[idx],noise[idx],
				   seed[idx]);
  } 
  return EXIT_SUCCESS;
}
int yoga_sensors::sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size, 
				 float *noise)
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx],ypos[idx],lambda[idx],mag[idx],size[idx],noise[idx],
				   1234*idx);
  } 
  return EXIT_SUCCESS;
}
int yoga_sensors::sensors_initgs(float *xpos,float *ypos,float *lambda, float *mag, long *size)
{
  for (size_t idx = 0; idx < (this->d_wfs).size(); idx++) {
    (this->d_wfs)[idx]->wfs_initgs(xpos[idx],ypos[idx],lambda[idx],mag[idx],size[idx],-1,1234);
  } 
  return EXIT_SUCCESS;
}
