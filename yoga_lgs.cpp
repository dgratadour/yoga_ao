#include <yoga_ao_utils.h>
#include <yoga_lgs.h>


yoga_lgs::yoga_lgs(long nvalid, long npix)
{
  this->nvalid  = nvalid;
  this->npix    = npix;
  this->nprof   = nprof;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  long *dims_data2 = new long[3];
  dims_data2[0] = 2; 
  long *dims_data3 = new long[4];
  dims_data3[0] = 3;

  dims_data1[1] = npix; 
  this->d_beam = new yoga_obj<float>(dims_data1);
  this->d_ftbeam = new yoga_obj<cuFloatComplex>(dims_data1);

  dims_data1[1] = nvalid; 
  this->d_azimuth = new yoga_obj<float>(dims_data1);

  dims_data2[1] = npix; dims_data2[2] = nvalid;
  this->d_prof2d = new yoga_obj<cuFloatComplex>(dims_data2);

  int mdims[1];
  mdims[0] = (int)dims_data2[1];
  cufftSafeCall(cufftPlanMany(&(this->d_prof2d->plan), 1 ,mdims,NULL,1,0,NULL,1,0,
			      CUFFT_C2C ,(int)dims_data2[2]));
  this->d_prof2d->fft_on = true;
  
  dims_data3[1] = npix; dims_data3[2] = npix; dims_data3[3] = nvalid;  
  this->d_lgskern   = new yoga_obj<float>(dims_data3);
  this->d_ftlgskern = new yoga_obj<cuFloatComplex>(dims_data3);

  int mdims2[2];
  mdims[0] = (int)dims_data3[1];
  mdims[1] = (int)dims_data3[2];
  cufftSafeCall( cufftPlanMany(&(this->d_ftlgskern->plan), 2 ,mdims,NULL,1,0,NULL,1,0,CUFFT_C2C ,
			       (int)dims_data3[3]));
  this->d_ftlgskern->fft_on = true;

  cudaExtent volumeSize = make_cudaExtent(npix,npix,nvalid);

  this->channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  cutilSafeCall( cudaMalloc3DArray(&(this->d_spotarray), &(this->channelDesc), volumeSize) );

  // prepare 3d cppy
  this->copyParams.srcPtr = make_cudaPitchedPtr((void*)(this->d_lgskern->d_data), volumeSize.width*sizeof(float), 
						volumeSize.width, volumeSize.height);
  copyParams.dstArray = this->d_spotarray;
  copyParams.extent   = volumeSize;
  copyParams.kind     = cudaMemcpyDeviceToDevice;
}

yoga_lgs::~yoga_lgs()
{
  delete this->d_prof1d;
  delete this->d_profcum;
  delete this->d_prof2d;
  delete this->d_ftbeam;
  delete this->d_beam;
  delete this->d_lgskern;
  delete this->d_ftlgskern;
  delete this->d_azimuth;
  delete this->d_doffaxis;
  cutilSafeCall(cudaFreeArray(this->d_spotarray));
}

int yoga_lgs::lgs_init(int nprof, float hg, float h0, float deltah, float pixsize, float *doffaxis, float *prof1d,
		       float *profcum,float *beam,cuFloatComplex *ftbeam,float *azimuth)
{
  this->nprof = nprof;
  this->hg = hg;
  this->h0 = h0;
  this->deltah = deltah;
  this->pixsize = pixsize;

  long *dims_data1 = new long[2];
  dims_data1[0] = 1;
  dims_data1[1] = this->nvalid;
  this->d_doffaxis = new yoga_obj<float>(dims_data1);
  this->d_doffaxis->host2device(doffaxis);

  dims_data1[1] = nprof+1; 
  this->d_prof1d = new yoga_obj<float>(dims_data1);
  this->d_profcum = new yoga_obj<float>(dims_data1);

  this->d_prof1d->host2device(prof1d);
  this->d_profcum->host2device(profcum);
  this->d_beam->host2device(beam);
  this->d_ftbeam->host2device(ftbeam);
  this->d_azimuth->host2device(azimuth);

  return EXIT_SUCCESS;
}

int yoga_lgs::load_prof(float *prof1d,float *profcum, float hg, float h0, float deltah)
{
  this->d_prof1d->host2device(prof1d);
  this->d_profcum->host2device(profcum);
  this->hg = hg;
  this->h0 = h0;
  this->deltah = deltah;

  return EXIT_SUCCESS;
}

int yoga_lgs::lgs_update(int device)
{
  interp_prof(this->d_prof2d->d_data,this->d_prof1d->d_data,this->d_profcum->d_data,
	       this->npix,this->d_doffaxis->d_data,this->hg,this->pixsize,this->h0,
	       this->deltah,this->nprof,this->d_prof2d->nb_elem,device);

  // convolution by beam
  // do fft on prof2d
  yoga_fft(this->d_prof2d->d_data,this->d_prof2d->d_data,1,this->d_prof2d->plan);

  // mult by beamft
  times_ftbeam(this->d_prof2d->d_data,this->d_ftbeam->d_data,this->npix, 
	       this->d_prof2d->nb_elem,device);

  // fft back
  yoga_fft(this->d_prof2d->d_data,this->d_prof2d->d_data,-1,this->d_prof2d->plan);

  // build final image
  // get abs of real and roll
  rollbeamexp(this->d_lgskern->d_data,this->d_prof2d->d_data,this->d_beam->d_data,
	  this->npix,this->d_lgskern->nb_elem,device);

  // rotate image and fill kernels ft
  /*
  lgs_rotate(this->d_ftlgskern->d_data,this->d_lgskern->d_data,this->npix,this->npix,
	     this->d_azimuth->d_data,0.0f,this->npix*this->npix*this->nvalid,device);
  */
  //same with textures
  rotate3d(this->d_ftlgskern->d_data,this->copyParams,this->d_spotarray,this->channelDesc,
	   this->npix,this->npix,this->d_azimuth->d_data,0.0f,this->npix*this->npix*this->nvalid,device);
  //prepare for wfs code
  yoga_fft(this->d_ftlgskern->d_data,this->d_ftlgskern->d_data,1,this->d_ftlgskern->plan);
  
  return EXIT_SUCCESS;
}

int yoga_lgs::load_kernels(float *lgskern,int device)
{
  this->d_lgskern->host2device(lgskern);
  cfillrealp(this->d_ftlgskern->d_data,this->d_lgskern->d_data,this->d_ftlgskern->nb_elem,device);
  yoga_fft(this->d_ftlgskern->d_data,this->d_ftlgskern->d_data,1,this->d_ftlgskern->plan);

  return EXIT_SUCCESS;
}

