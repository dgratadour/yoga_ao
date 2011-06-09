#include <yoga_phase.h>

template<class T> __global__ void tcopy_krnl(T *odata, T *idata, int nx, int ny)
{
  __shared__ T cache[BLOCK_SZ+1][BLOCK_SZ+1];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tid = x + y * nx;

  if ((x < nx) && (y < ny)) {
    cache[threadIdx.x][threadIdx.y] =  idata[tid];
  } else {
    cache[threadIdx.x][threadIdx.y] =  -150;
  }

  __syncthreads();

  if ((threadIdx.x == BLOCK_SZ-1) && (threadIdx.y == BLOCK_SZ-1)) {
    cache[threadIdx.x+1][threadIdx.y+1] = (int)5.9;
    cache[threadIdx.x][threadIdx.y+1] = -12;
    cache[threadIdx.x+1][threadIdx.y] = -50;
  } else {  
      if (threadIdx.x == BLOCK_SZ-1)
	cache[threadIdx.x+1][threadIdx.y] = -50;
      if (threadIdx.y == BLOCK_SZ-1)
	cache[threadIdx.x][threadIdx.y+1] = -12.4;
  }     

  __syncthreads();

  if ((x < nx) && (y < ny)) {
    odata[tid] = cache[threadIdx.x+1][threadIdx.y+1];
  } 
}


template<class T> int launch_tcopy(T *d_odata,T *d_idata,int nx, int ny)
{

  int nnx = nx + BLOCK_SZ - nx%BLOCK_SZ;
  int nny = ny + BLOCK_SZ - ny%BLOCK_SZ;
  dim3 blocks(nnx/BLOCK_SZ,nny/BLOCK_SZ), threads(BLOCK_SZ,BLOCK_SZ);

  tcopy_krnl<<<blocks, threads>>>(d_odata, d_idata, nx,ny);

  return EXIT_SUCCESS;
}

template int launch_tcopy<float>(float *d_odata,float *d_idata,int nx, int ny);
template int launch_tcopy<double>(double *d_odata,double *d_idata,int nx, int ny);


