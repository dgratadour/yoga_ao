#include <yoga_wfs.h>

__global__ void krnl_fillcamplipup(cuFloatComplex *amplipup, float *phase,float *offset, float *mask, int *indx, int Nfft, int Npup, int npup, int N)
{
  int nim,idim,idimx,idimy,idx;
  int tid   = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim   = tid / Npup;
    idim  = tid - nim * Npup;
    idimx = idim % npup;
    idimy = idim / npup;
    idx   = idimx + idimy * Nfft + nim * Nfft * Nfft;
    amplipup[idx].x = (cosf(phase[indx[tid]]-offset[idim]))*mask[indx[tid]];
    amplipup[idx].y = (sinf(phase[indx[tid]]-offset[idim]))*mask[indx[tid]];
    //amplipup[idx].x = mask[indx[tid]] * cosf(-offset[idim]);
    //amplipup[idx].y = mask[indx[tid]] * sinf(-offset[idim]);
    tid  += blockDim.x * gridDim.x;
  }
}

int fillcamplipup(cuFloatComplex *amplipup, float *phase, float *offset, float *mask, int *indx, int Nfft, int Npup, int Nsub, int npup, int device)
// here amplipup is a cube of data of size nfft x nfft x nsubap
// phase is an array of size pupdiam x pupdiam
// offset is an array of size pdiam x pdiam
// mask is an array of size pupdiam x pupdiam
// indx is an array of size pdiam x pdiam x nsubap
// number of thread required : pdiam x pdiam x nsubap
// Npup = pdiam x pdiam
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  
  int N = Npup * Nsub;

  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_fillcamplipup<<<grid, threads>>>(amplipup,phase,offset,mask,indx,Nfft,Npup,npup,N);
  cutilCheckMsg("fillcamplipup_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void krnl_fillbincube(float *bcube, float *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin, int N)
{
  /*
    indx is an array nrebin^2 * npix^2
    it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
    Npix = npix x npix
   */
  int npix,nsubap,nrebin;
  int tid     = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nsubap  = tid / Npix;
    npix    = tid % Npix;
    for (int i=0;i<Nrebin;i++) {
      nrebin = indxpix[i+npix*Nrebin];
      bcube[tid] += hrimage[nrebin + Nfft*nsubap];
    }
    tid  += blockDim.x * gridDim.x;
  }
}

int fillbincube(float *bcube, float *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin, int Nsub, int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  
  int N = Npix * Nsub;
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_fillbincube<<<grid, threads>>>(bcube,hrimage,indxpix,Nfft,Npix,Nrebin,N);
  cutilCheckMsg("fillbincube_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void krnl_fillbinimg(float *bimage, float *bcube, int *indxbin, int N)
{
  /*
    indx is an array nrebin^2 * npix^2
    it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
    Npix = npix x npix
   */
  int tid     = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    bimage[indxbin[tid]] = bcube[tid];
    tid  += blockDim.x * gridDim.x;
  }
}

int fillbinimg(float *bimage, float *bcube, int *indxbin, int Npix, int Nsub, int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  
  int N = Npix * Nsub;
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_fillbinimg<<<grid, threads>>>(bimage,bcube,indxbin,N);
  cutilCheckMsg("fillbinimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

template<class T> __global__ void krnl_indexfill(T *odata, T *idata, int *indx,int ntot, int Ntot, int N)
{
  int nim, idim;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / ntot;
    idim = tid - (nim * ntot);
    odata[indx[idim] + (nim * Ntot)] = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int indexfill(float *d_odata,float *d_idata,int *indx,int nx, int Nx, int N,int device)
{

  int ntot = nx * nx;
  int Ntot = Nx * Nx;
  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
    
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  krnl_indexfill<<<grid, threads>>>(d_odata, d_idata,indx, ntot, Ntot, N);

  return EXIT_SUCCESS;
}

__global__ void cfillrealp(cuFloatComplex *odata, float *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid].x = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int launch_cfillrealp(cuFloatComplex *d_odata,float *d_idata,int N,int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
    
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  cfillrealp<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}

__global__ void conv_krnl(cuFloatComplex *odata,cuFloatComplex *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex tmp;

  while (tid < N) {
    tmp.x = idata[tid].x*odata[tid].x+idata[tid].y*odata[tid].y;
    tmp.y = -1.0f*idata[tid].y*odata[tid].x+idata[tid].x*odata[tid].y;
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

int launch_conv_krnl(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
    
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  conv_krnl<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}

__global__ void cgetrealp(float *odata, cuFloatComplex *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[tid].x;
    tid += blockDim.x * gridDim.x;
  }
}

int launch_cgetrealp(float *d_odata,cuFloatComplex *d_idata,int N,int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
    
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (N + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);
  //  dim3 grid(128), threads(128);

  cgetrealp<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}

