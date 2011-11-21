#include <yoga_ao_utils.h>

__global__ void cfillrealp_krnl(cuFloatComplex *odata, float *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid].x = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int cfillrealp(cuFloatComplex *d_odata,float *d_idata,int N,int device)
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

  cfillrealp_krnl<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}


__global__ void cgetrealp_krnl(float *odata, cuFloatComplex *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[tid].x;
    tid += blockDim.x * gridDim.x;
  }
}

int cgetrealp(float *d_odata,cuFloatComplex *d_idata,int N,int device)
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

  cgetrealp_krnl<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}


template<class T> __global__ void roll_krnl(T *idata,int N,int  M,int Ntot)
{

  int tidt = threadIdx.x + blockIdx.x * blockDim.x;
  int nim  = tidt / Ntot;

  int tid = tidt - nim * Ntot;

  while (tid < Ntot) {
   
    
   int x=tid%N;
   int y=tid/N;

   int  xx=(x+N/2)%N;
   int  yy=(y+M/2)%M;
   int  tid2=xx+yy*N;

   __shared__ T tmp;
   tmp=idata[tid + nim * (N*M)];
   idata[tid + nim * (N*M)]=idata[tid2 + nim * (N*M)];
   idata[tid2 + nim * (N*M)]=tmp;

    tid += blockDim.x * gridDim.x;
  }
}

template<class T> int roll(T *idata,int N,int M,int nim)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, 0);
    
  long Ntot = N * M * nim;
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*8;
  int nThreads = (Ntot/2 + nBlocks -1)/nBlocks;

  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (Ntot/2 + nThreads -1)/nThreads;
  }

  dim3 grid(nBlocks), threads(nThreads);

  roll_krnl<<<grid, threads>>>(idata,N,M,Ntot/2);

  return EXIT_SUCCESS;

}

template int roll<float>(float *idata,int N,int M,int nim);

template int roll<double>(double *idata,int N,int M,int nim);

template int roll<cuFloatComplex>(cuFloatComplex *idata,int N,int M,int nim);


__global__ void abs2_krnl(float *odata, cuFloatComplex *idata, int N)
{
  __shared__ cuFloatComplex cache;

  int tid   = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < N) {
    cache = idata[tid];
    odata[tid] = cache.x *cache.x +cache.y *cache.y;
  } 
}


int abs2(float *d_odata, cuFloatComplex *d_idata, int N, int device)
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

  abs2_krnl<<<grid, threads>>>(d_odata,d_idata,N);
  cutilCheckMsg("abs2_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void abs2c_krnl(cuFloatComplex *odata, cuFloatComplex *idata, int N)
{
  __shared__ cuFloatComplex cache;
  int tid   = threadIdx.x + blockIdx.x * blockDim.x;


  if (tid < N) {
    cache = idata[tid];
    odata[tid].x = cache.x *cache.x +cache.y *cache.y;
    odata[tid].y = 0.0;
  } 
}


int abs2c(cuFloatComplex *d_odata, cuFloatComplex *d_idata, int N, int device)
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

  abs2c_krnl<<<grid, threads>>>(d_odata,d_idata,N);
  cutilCheckMsg("abs2c_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}


__global__ void subapnorm_krnl(float *odata, float *idata, float *fact, float *norm, float nphot,int n, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[tid] * fact[tid/n] / norm[tid/n] * nphot;
    tid += blockDim.x * gridDim.x;
  }
}

int subap_norm(float *d_odata,float *d_idata,float *fact,float *norm,float nphot,int n, int N,int device)
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

  subapnorm_krnl<<<grid, threads>>>(d_odata, d_idata, fact, norm, nphot, n, N);

   return EXIT_SUCCESS;
}


__global__ void krnl_fillindx(float *odata, float *idata, int *indx, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] = idata[indx[tid]];
    tid += blockDim.x * gridDim.x;
  }
}

int fillindx(float *d_odata,float *d_idata,int *indx,int N, int device)
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

  krnl_fillindx<<<grid, threads>>>(d_odata, d_idata,indx, N);

  return EXIT_SUCCESS;
}

__global__ void fillarr2d_krnl(float *odata, float *idata, int tidx0, int Ncol,int NC, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

  while (tid < N) {
    if (Ncol > 1) tidB = tidx0 + (tid/Ncol)*NC + (tid%Ncol);
    else tidB = tidx0 + tid*NC;
    odata[tidB] = idata[tid];
    tid += blockDim.x * gridDim.x;
  }
}

int fillarr2d(float *d_odata,float *d_idata,int x0, int Ncol,int NC, int N,int device)
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

  fillarr2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol,NC,N);

  cutilCheckMsg("fillarr2d_kernel<<<>>> execution failed\n");

   return EXIT_SUCCESS;
}

__global__ void getarr2d_krnl(float *odata, float *idata, int tidx0, int Ncol,int NC, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tidB;

  while (tid < N) {
    if (Ncol > 1) tidB = tidx0 + (tid/Ncol)*NC + (tid%Ncol);
    else tidB = tidx0 + tid*NC;
    odata[tid] = idata[tidB];
    tid += blockDim.x * gridDim.x;
  }
}

int getarr2d(float *d_odata,float *d_idata,int x0, int Ncol,int NC, int N, int device)
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

  getarr2d_krnl<<<grid, threads>>>(d_odata, d_idata, x0, Ncol,NC,N);

  cutilCheckMsg("getarr2d_kernel<<<>>> execution failed\n");

   return EXIT_SUCCESS;
}

__global__ void addai_krnl(float *odata, float* idata, int i, int sgn, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    if (sgn == 1) odata[tid] += idata[i];
    else odata[tid] -= idata[i];
    tid += blockDim.x * gridDim.x;
  }
}

int addai(float *d_odata,float *i_data, int i,int sgn, int N, int device)
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

  addai_krnl<<<grid, threads>>>(d_odata, i_data, i, sgn, N);

  cutilCheckMsg("plusai_kernel<<<>>> execution failed\n");

   return EXIT_SUCCESS;
}


