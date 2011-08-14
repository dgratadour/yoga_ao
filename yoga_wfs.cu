#include <yoga_wfs.h>
bool isPow2(unsigned int x)
{
    return ((x&(x-1))==0);
}

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T> struct SharedMemory
{
    __device__ inline operator       T*()
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }

    __device__ inline operator const T*() const
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }
};

// specialize for double to avoid unaligned memory 
// access compile errors
template<> struct SharedMemory<double>
{
    __device__ inline operator       double*()
    {
        extern __shared__ double __smem_d[];
        return (double*)__smem_d;
    }

    __device__ inline operator const double*() const
    {
        extern __shared__ double __smem_d[];
        return (double*)__smem_d;
    }
};



__global__ void fillcamplipup_krnl(cuFloatComplex *amplipup, float *phase,float *offset, float *mask, int *indx, int Nfft, int Npup, int npup, int N)
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

  fillcamplipup_krnl<<<grid, threads>>>(amplipup,phase,offset,mask,indx,Nfft,Npup,npup,N);
  cutilCheckMsg("fillcamplipup_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}



__global__ void fillbincube_krnl(float *bcube, float *hrimage, int *indxpix, int Nfft, int Npix, int Nrebin, int N)
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

  fillbincube_krnl<<<grid, threads>>>(bcube,hrimage,indxpix,Nfft,Npix,Nrebin,N);
  cutilCheckMsg("fillbincube_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}



__global__ void fillbinimg_krnl(float *bimage, float *bcube, int *indxbin, int N)
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

__global__ void addbinimg_krnl(float *bimage, float *bcube, int *indxbin, int N)
{
  /*
    indx is an array nrebin^2 * npix^2
    it gives the nrebin x nrebin pixels in the hrimage per npix x npix pixels of the subap
    Npix = npix x npix
   */
  int tid     = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    bimage[indxbin[tid]] += bcube[tid];
    tid  += blockDim.x * gridDim.x;
  }
}

int fillbinimg(float *bimage, float *bcube, int *indxbin, int Npix, int Nsub, bool add, int device)
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

  if (add)
    addbinimg_krnl<<<grid, threads>>>(bimage,bcube,indxbin,N);
  else 
    fillbinimg_krnl<<<grid, threads>>>(bimage,bcube,indxbin,N);

  cutilCheckMsg("fillbinimg_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}



template<class T> __global__ void indexfill_krnl(T *odata, T *idata, int *indx,int ntot, int Ntot, int N)
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

  indexfill_krnl<<<grid, threads>>>(d_odata, d_idata,indx, ntot, Ntot, N);

  return EXIT_SUCCESS;
}



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



__global__ void conv_krnl(cuFloatComplex *odata,cuFloatComplex *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  cuFloatComplex tmp;

  while (tid < N) {
    tmp.x = idata[tid].x*odata[tid].x-idata[tid].y*odata[tid].y;
    tmp.y = idata[tid].y*odata[tid].x+idata[tid].x*odata[tid].y;
    odata[tid] = tmp;
    tid += blockDim.x * gridDim.x;
  }
}

int convolve(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device)
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



/*
    This version uses sequential addressing -- no divergence or bank conflicts.
*/
template <class T> __device__ void reduce_krnl(T *sdata, int size, int n)
{
  if (!((size&(size-1))==0)) {
    unsigned int s;
    if (size %2 != 0) s = size/2+1;
    else s = size/2;
    unsigned int s_old = size;
    while (s>0) {
      if ((n < s) && (n + s < s_old)) {
	sdata[n] += sdata[n + s];
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2*s < s_old) && (s!=0)) s += 1;
    }
  } else {
    // do reduction in shared mem
    for(unsigned int s=size/2; s>0; s>>=1) {	
      if (n < s) {
	sdata[n] += sdata[n + s];
      }
      __syncthreads();
    }
  }
}


template <class T> __device__ void scanmax_krnl(T *sdata, int *values,int size, int n)
{
  if (!((size&(size-1))==0)) {
    unsigned int s;
    if (size %2 != 0) s = size/2+1;
    else s = size/2;
    unsigned int s_old = size;
    while (s>0) {
      if ((n < s) && (n + s < s_old)) {
	if (sdata[n] < sdata[n + s]) {
	  values[n] = n + s;
	  sdata[n] = sdata[n+s];
	}
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2*s < s_old) && (s!=0)) s += 1;
    }
  } else {
    // do reduction in shared mem
    for(unsigned int s=size/2; s>0; s>>=1) {	
      if (n < s) {
	if (sdata[n] < sdata[n + s]) {
	  values[n] = n + s;
	  sdata[n] = sdata[n+s];
	}
      }
      __syncthreads();
    }
  }
}


template <class T> 
__device__ inline void mswap(T & a, T & b)
{
    T tmp = a;
    a = b;
    b = tmp;
}


template <class T> __device__ inline void sortmax_krnl(T *sdata, unsigned int *values,int size, int n)
{

  if (!((size&(size-1))==0)) {
    unsigned int s;
    if (size %2 != 0) s = size/2+1;
    else s = size/2;
    unsigned int s_old = size;
    while (s>0) {
      if ((n < s) && (n + s < s_old)) {
	if (sdata[n] < sdata[n + s]) {
	  mswap(values[n],values[n+s]);
	  mswap(sdata[n],sdata[n+s]);
	}
      }
      __syncthreads();
      s_old = s;
      s /= 2;
      if ((2*s < s_old) && (s!=0)) s += 1;
    }
  } else {
    // do reduction in shared mem
    for(unsigned int s=size/2; s>0; s>>=1) {	
      if (n < s) {
	if (sdata[n] < sdata[n + s]) {
	  mswap(values[n],values[n+s]);
	  mswap(sdata[n],sdata[n+s]);
	}
      }
      __syncthreads();
    }
  }
}

template <class T> __global__ void reduce2(T *g_idata, T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata[tid] = (i < n) ? g_idata[i] : 0;
    
    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class T> __global__ void reduce2(T *g_idata, T *g_odata, T thresh, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    if (i < n) sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] : 0;
    
    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);
    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

template <class T> __global__ void sortmax(T *g_idata, T *g_odata, int *values, int nmax)
{
  extern __shared__ uint svalues[];
  T *sdata = (T*)&svalues[blockDim.x];
  
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
  
  svalues[tid] = tid;
  sdata[tid]   = g_idata[i];
  
  __syncthreads();
  
  for (int cc=0;cc<nmax;cc++) {
    
    if (tid >= cc) sortmax_krnl(&(sdata[cc]),&(svalues[cc]),blockDim.x-cc,tid-cc);
    
    __syncthreads();
    
  } 
  
  if (tid < nmax) {
    g_odata[nmax*blockIdx.x+tid] = sdata[tid];
    values[nmax*blockIdx.x+tid]  = svalues[tid];
  }
  __syncthreads();
 
}

template <class T> __global__ void centroid_max(T *g_idata, T *g_odata, int n, int nmax, int nsub)
{
  extern __shared__ uint svalues[];
  T *sdata = (T*)&svalues[blockDim.x];
  T subsum;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
  
  svalues[tid] = tid;
  sdata[tid]   = g_idata[i];
  
  __syncthreads();
  
  for (int cc=0;cc<nmax;cc++) {
    
    if (tid >= cc) sortmax_krnl(&(sdata[cc]),&(svalues[cc]),blockDim.x-cc,tid-cc);
    
    __syncthreads();
    
  } 
  // at this point the nmax first elements of sdata are the nmax brightest
  // pixels and the nmax first elements of svalues are their positions

  // first copy the brightest values out for reduction  
  if ((tid >= nmax) && (tid < 2*nmax)) sdata[tid] = sdata[tid-nmax];

  __syncthreads();

  reduce_krnl(sdata,nmax,tid);

  // get the sum per subap  
  if (tid == 0) subsum = sdata[tid];

  // put back the brightest pixels values 
  if ((tid >= nmax) && (tid < 2*nmax)) sdata[tid-nmax] = sdata[tid];

  __syncthreads();

  // compute the centroid on the first part of the array
  if (tid < nmax) sdata[tid] *= ((svalues[tid] % n) + 1);
  // x centroid
  __syncthreads();
  reduce_krnl(sdata,nmax,tid);

  if (tid == 0) g_odata[blockIdx.x] = sdata[tid]/subsum;

  // put back the brightest pixels values 
  if ((tid >= nmax) && (tid < 2*nmax)) sdata[tid-nmax] = sdata[tid];

  __syncthreads();

  // compute the centroid on the first part of the array
  if (tid < nmax) sdata[tid] *= (svalues[tid] / n + 1);
  // y centroid
  __syncthreads();
  reduce_krnl(sdata,nmax,tid);

  if (tid == 0) g_odata[blockIdx.x+nsub] = sdata[tid]/subsum;
 
}


template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
    reduce2<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size);

    cutilCheckMsg("reduce_kernel<<<>>> execution failed\n");
}
template void subap_reduce<float>(int size, int threads, int blocks, float *d_idata, float *d_odata);

template void subap_reduce<double>(int size, int threads, int blocks, double *d_idata, double *d_odata);


template <class T> void subap_reduce(int size, int threads, int blocks, T *d_idata, T *d_odata, T thresh)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
    reduce2<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, thresh, size);

    cutilCheckMsg("reduce_kernel<<<>>> execution failed\n");
}
template void subap_reduce<float>(int size, int threads, int blocks, float *d_idata, float *d_odata, float thresh);

template void subap_reduce<double>(int size, int threads, int blocks, double *d_idata, double *d_odata, double thresh);




template <class T> void subap_sortmax(int size, int threads, int blocks, T *d_idata, T *d_odata, int *values, int nmax)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * (sizeof(T) + sizeof(uint)) : threads * (sizeof(T) + sizeof(uint));
    sortmax<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, values, nmax);

    cutilCheckMsg("sortmax_kernel<<<>>> execution failed\n");
}
template void subap_sortmax<float>(int size, int threads, int blocks, float *d_idata, float *d_odata, int *values, int nmax);

template void subap_sortmax<double>(int size, int threads, int blocks, double *d_idata, double *d_odata, int *values, int nmax);


template <class T> void subap_centromax(int threads, int blocks, T *d_idata, T *d_odata, int npix, int nmax)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * (sizeof(T) + sizeof(uint)) : threads * (sizeof(T) + sizeof(uint));
    centroid_max<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata,npix, nmax,blocks);

    cutilCheckMsg("sortmax_kernel<<<>>> execution failed\n");
}
template void subap_centromax<float>( int threads, int blocks, float *d_idata, float *d_odata, int npix, int nmax);

template void subap_centromax<double>(int threads, int blocks, double *d_idata, double *d_odata, int npix, int nmax);




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




template <class T> __global__ void reduce_phasex(T *g_idata, T *g_odata, T *offset, int *indx, unsigned int n, T alpha)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i   = blockIdx.x * n * n;
        
    sdata[tid] = (g_idata[indx[i+tid*n+n-1]]-offset[tid*n+n-1])-
      (g_idata[indx[i+tid*n]]-offset[tid*n]);
    
    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) {
      g_odata[blockIdx.x] = sdata[0] / n * alpha;
    }
}

template <class T> __global__ void reduce_phasey(T *g_idata, T *g_odata, T *offset, int *indx, unsigned int n, T alpha)
// FIXME
// full parallelization would require to use 4 x threads
// instead of doing 4 operations per threads
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * n * n + threadIdx.x;
    
      sdata[tid] = (g_idata[indx[i+(n-1)*n]]-offset[(n-1)*n+tid])-
	(g_idata[indx[i]]-offset[tid]);

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) {
      g_odata[blockIdx.x] = sdata[0]/n*alpha;
    }
}

template <class T> void phase_reduce(int threads, int blocks, T *d_idata, T *d_odata, int *indx, T *offset, T alpha)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
    reduce_phasex<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, offset, indx, threads, alpha);

    cutilCheckMsg("reduce_phasex_kernel<<<>>> execution failed\n");

    reduce_phasey<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, &(d_odata[blocks]), offset, indx, threads, alpha);

    cutilCheckMsg("reduce_phasey_kernel<<<>>> execution failed\n");
}

template void phase_reduce<float>(int threads, int blocks, float *d_idata, float *d_odata, int *indx, float *offset, float alpha);

template void phase_reduce<double>(int threads, int blocks, double *d_idata, double *d_odata, int *indx, double *offset, double alpha);




template <class T> __global__ void derive_phasex(T *g_idata, T *g_odata, int *indx, T *offset, T *mask, T alpha, unsigned int n, unsigned int N, float *fluxPerSub)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    if (i < N) {
      if (tid % n == 0) {//start of a column
	sdata[tid] = ((g_idata[indx[i+1]]-offset[tid+1]) - 
		      (g_idata[indx[i]]-offset[tid])) * mask[indx[i]];
      } else {
	if ((tid+1) % n == 0) {//end of a column
	  sdata[tid] = ((g_idata[indx[i]]-offset[tid]) -
			(g_idata[indx[i-1]]-offset[tid-1])) * mask[indx[i]];
	} else 
	  sdata[tid] = ((g_idata[indx[i+1]]-offset[tid+1]) - 
			(g_idata[indx[i-1]]-offset[tid-1])) * mask[indx[i]];
      }
    }

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0] / n / 2 * alpha / fluxPerSub[blockIdx.x];
}

template <class T> __global__ void derive_phasey(T *g_idata, T *g_odata, int *indx, T *offset, T *mask, T alpha, unsigned int n, unsigned int N, float *fluxPerSub)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    if (i < N) {
      if (tid < n) {//start of a column
	sdata[tid] = ((g_idata[indx[i+n]]-offset[tid+n]) - 
		      (g_idata[indx[i]]-offset[tid])) * mask[indx[i]];
      } else {
	if (tid >= n*(n-1)) {//end of a column
	  sdata[tid] = ((g_idata[indx[i]]-offset[tid]) -
			(g_idata[indx[i-n]]-offset[tid-n])) * mask[indx[i]];
	} else 
	  sdata[tid] = ((g_idata[indx[i+n]]-offset[tid+n]) - 
			(g_idata[indx[i-n]]-offset[tid-n])) * mask[indx[i]];
      }
    }

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0] / n / 2 * alpha / fluxPerSub[blockIdx.x];
}

template <class T> void phase_derive(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, int *indx, T *offset, T *mask, T alpha, float *fluxPerSub)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
    derive_phasex<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, indx, offset, mask, alpha, n, size,fluxPerSub);

    cutilCheckMsg("phase_derivex_kernel<<<>>> execution failed\n");

    derive_phasey<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, &(d_odata[blocks]), indx, offset, mask, alpha, n, size,fluxPerSub);

    cutilCheckMsg("phase_derivey_kernel<<<>>> execution failed\n");
 }

template void phase_derive<float>(int size, int threads, int blocks, int n, float *d_idata, float *d_odata, int *indx, float *offset, float *mask,float alpha, float *fluxPerSub);

template void phase_derive<double>(int size, int threads, int blocks, int n, double *d_idata, double *d_odata, int *indx, double *offset, double *mask,double alpha, float *fluxPerSub);



template <class T> __global__ void centroidx(T *g_idata, T *g_odata, T *alpha, unsigned int n, unsigned int N)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata[tid] = (i < N) ? g_idata[i] * ((tid % n) +1) : 0;

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0]*1.0/(alpha[blockIdx.x]+1.e-6);
}

template <class T> __global__ void centroidy(T *g_idata, T *g_odata, T *alpha, unsigned int n, unsigned int N)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    sdata[tid] = (i < N) ? g_idata[i] * ((tid / n) +1) : 0;

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0]*1.0/(alpha[blockIdx.x]+1.e-6);
}

template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, T *alpha)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
    centroidx<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, alpha, n, size);

    cutilCheckMsg("centroidx_kernel<<<>>> execution failed\n");

    centroidy<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, &(d_odata[blocks]), alpha, n, size);

    cutilCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int size, int threads, int blocks, int n,float *d_idata, float *d_odata,float *alpha);

template void get_centroids<double>(int size, int threads, int blocks, int n,double *d_idata, double *d_odata, double *alpha);


template <class T> __global__ void centroidx(T *g_idata, T *g_odata, T *alpha, T thresh, unsigned int n, unsigned int N)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    if (i < N) sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] * ((tid % n) +1) : 0;

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0]*1.0/(alpha[blockIdx.x]+1.e-6);
}

template <class T> __global__ void centroidy(T *g_idata, T *g_odata, T *alpha, T thresh, unsigned int n, unsigned int N)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
    
    if (i < N) sdata[tid] = (g_idata[i] > thresh) ? g_idata[i] * ((tid / n) +1) : 0;

    __syncthreads();

    reduce_krnl(sdata,blockDim.x,tid);

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0]*1.0/(alpha[blockIdx.x]+1.e-6);
}

template <class T> void get_centroids(int size, int threads, int blocks, int n, T *d_idata, T *d_odata, T *alpha, T thresh)
{
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);

    // when there is only one warp per block, we need to allocate two warps 
    // worth of shared memory so that we don't index shared memory out of bounds
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);
    centroidx<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, alpha, thresh, n, size);

    cutilCheckMsg("centroidx_kernel<<<>>> execution failed\n");

    centroidy<T><<< dimGrid, dimBlock, smemSize >>>(d_idata, &(d_odata[blocks]), alpha, thresh, n, size);

    cutilCheckMsg("centroidy_kernel<<<>>> execution failed\n");
}

template void get_centroids<float>(int size, int threads, int blocks, int n,float *d_idata, float *d_odata,float *alpha, float thresh);

template void get_centroids<double>(int size, int threads, int blocks, int n,double *d_idata, double *d_odata, double *alpha, double thresh);


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


__global__ void fillcorr_krnl(cuFloatComplex *d_out, float *d_in,int npix_in,int npix_out,int N)
{
  int nim,npix,idx;
  int tid   = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / npix_in;
    npix = tid % npix_in;
    idx = nim * npix_out + npix;
    d_out[idx].x = d_in[tid];
    d_out[idx].y = 0.0;
    tid  += blockDim.x * gridDim.x;
  }
}

int fill_corr(cuFloatComplex *d_out, float *d_in, int npix_in, int npix_out, int N, int device)
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

  fillcorr_krnl<<<grid, threads>>>(d_out,d_in,npix_in,npix_out,N);
  cutilCheckMsg("fillcorr_kernel<<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void fillval_krnl(cuFloatComplex *d_out, float val,int npix_in,int npix_out,int N)
{
  int nim,npix,idx;
  int tid   = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    nim = tid / npix_in;
    npix = tid % npix_in;
    idx = nim * npix_out + npix;
    d_out[idx].x = val;
    d_out[idx].y = 0.0;
    tid  += blockDim.x * gridDim.x;
  }
}

int fillval_corr(cuFloatComplex *d_out, float val, int npix_in, int npix_out, int N, int device)
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

  fillval_krnl<<<grid, threads>>>(d_out,val,npix_in,npix_out,N);
  cutilCheckMsg("fillcorr_kernel<<<>>> execution failed\n");

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


__global__ void corr_krnl(cuFloatComplex *odata,cuFloatComplex *idata, int N)
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

int correl(cuFloatComplex *d_odata,cuFloatComplex *d_idata,int N,int device)
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

  corr_krnl<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}

__global__ void corrnorm_krnl(float *odata,float *idata, int N)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < N) {
    odata[tid] /= (idata[tid]+1.e-6);
    tid += blockDim.x * gridDim.x;
  }
}

int corr_norm(float *d_odata,float *d_idata,int N,int device)
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

  corrnorm_krnl<<<grid, threads>>>(d_odata, d_idata, N);

   return EXIT_SUCCESS;
}
