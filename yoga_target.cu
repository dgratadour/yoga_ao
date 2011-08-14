#include <yoga_target.h>

__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,float *xref, float *yref, int Nx)
{
  __shared__ float cache[BLOCK_SZ+1][BLOCK_SZ+1];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tido;

  int iref,jref,tidi;

  float xshift,yshift,wx1,wx2,wy1,wy2;

  if ((x < nx) && (y < ny)) {
    iref = (int)xref[x];
    jref = (int)yref[y];
    tidi = iref + jref * Nx;
    cache[threadIdx.x][threadIdx.y] =  idata[tidi];
  } 

  if ((x == nx-1) && (y < ny-1)) {
    cache[threadIdx.x+1][threadIdx.y] = idata[tidi+1];
    if (threadIdx.y == BLOCK_SZ-1) {
      cache[threadIdx.x][threadIdx.y+1] = idata[tidi+Nx];
      cache[threadIdx.x+1][threadIdx.y+1] = idata[tidi+Nx+1];
    }
  }

  if ((x < nx-1) && (y == ny-1)) {
    cache[threadIdx.x][threadIdx.y+1] = idata[tidi+Nx];
    if (threadIdx.x == BLOCK_SZ-1) {
      cache[threadIdx.x+1][threadIdx.y]   = idata[tidi+1];
      cache[threadIdx.x+1][threadIdx.y+1] = idata[tidi+Nx+1];
    }
  }

  if ((x == nx-1) && (y == ny-1)) {
    cache[threadIdx.x][threadIdx.y+1] = idata[tidi+Nx];
    cache[threadIdx.x+1][threadIdx.y] = idata[tidi+1];
    cache[threadIdx.x+1][threadIdx.y+1] = idata[tidi+Nx+1];
  }

  __syncthreads();
  
  if ((x < nx) && (y < ny)) {
    if ((threadIdx.x == BLOCK_SZ-1) && (threadIdx.y == BLOCK_SZ-1)) {
      cache[threadIdx.x+1][threadIdx.y+1] = idata[tidi+Nx+1];
      cache[threadIdx.x][threadIdx.y+1] = idata[tidi+Nx];
      cache[threadIdx.x+1][threadIdx.y] = idata[tidi+1];
    } else {  
      if (threadIdx.x == BLOCK_SZ-1)
	cache[threadIdx.x+1][threadIdx.y] = idata[tidi+1];
      if (threadIdx.y == BLOCK_SZ-1)
	cache[threadIdx.x][threadIdx.y+1] = idata[tidi+Nx];
    }     
  }

  __syncthreads();

  if ((x < nx) && (y < ny)) {
    tido = x + y * nx;

    xshift = xref[x] - iref;
    yshift = yref[y] - jref;

    wx1 = (1.0f - xshift);
    wx2 = xshift;
    wy1 = (1.0f - yshift);
    wy2 = yshift;

    odata[tido] += (wx1 * wy1 * cache[threadIdx.x][threadIdx.y] +
		    wx2 * wy1 * cache[threadIdx.x+1][threadIdx.y] +
		    wx1 * wy2 * cache[threadIdx.x][threadIdx.y+1] +
		    wx2 * wy2 * cache[threadIdx.x+1][threadIdx.y+1]);
  } 
}


int target_raytrace(float *d_odata,float *d_idata,int nx, int ny,float *xref, float *yref, int Nx)
{

  int nnx = nx + BLOCK_SZ - nx%BLOCK_SZ; // find next multiple of BLOCK_SZ
  int nny = ny + BLOCK_SZ - ny%BLOCK_SZ;
  dim3 blocks(nnx/BLOCK_SZ,nny/BLOCK_SZ), threads(BLOCK_SZ,BLOCK_SZ);

  raytrace_krnl<<<blocks, threads>>>(d_odata, d_idata, nx,ny,xref,yref,Nx);

  cutilCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void fillampli_krnl(cuFloatComplex *odata, float *idata, float *mask, int nx, int ny, int Nx)
{
  __shared__ cuFloatComplex cache[BLOCK_SZ][BLOCK_SZ];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tid = x + y * nx;

  if ((x < nx) && (y < ny)) {
    (cache[threadIdx.x][threadIdx.y]).x =  mask[tid]*cosf(idata[tid]);
    (cache[threadIdx.x][threadIdx.y]).y =  mask[tid]*sinf(idata[tid]);
  } else {
    (cache[threadIdx.x][threadIdx.y]).x =  0.0f;
    (cache[threadIdx.x][threadIdx.y]).y =  0.0f;
  }

  __syncthreads();

  tid = x + y * Nx;

  if ((x < nx) && (y < ny)) {
    odata[tid] = cache[threadIdx.x][threadIdx.y];
  } 
}


int fillampli(cuFloatComplex *d_odata,float *d_idata, float *mask,int nx, int ny, int Nx)
{

  int nnx = nx + BLOCK_SZ - nx%BLOCK_SZ;
  int nny = ny + BLOCK_SZ - ny%BLOCK_SZ;
  dim3 blocks(nnx/BLOCK_SZ,nny/BLOCK_SZ), threads(BLOCK_SZ,BLOCK_SZ);

  cutilCheckMsg("fillampli_kernel<<<>>> execution failed\n");
  fillampli_krnl<<<blocks, threads>>>(d_odata, d_idata, mask, nx,ny,Nx);

  return EXIT_SUCCESS;
}

__global__ void fillpupil_krnl(cuFloatComplex *odata, float *mask, int nx, int ny, int Nx)
{
  __shared__ cuFloatComplex cache[BLOCK_SZ][BLOCK_SZ];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tid = x + y * nx;

  if ((x < nx) && (y < ny)) {
    (cache[threadIdx.x][threadIdx.y]).x =  mask[tid];
    (cache[threadIdx.x][threadIdx.y]).y =  0.0f;
  } else {
    (cache[threadIdx.x][threadIdx.y]).x =  0.0f;
    (cache[threadIdx.x][threadIdx.y]).y =  0.0f;
  }

  __syncthreads();

  tid = x + y * Nx;

  if ((x < nx) && (y < ny)) {
    odata[tid] = cache[threadIdx.x][threadIdx.y];
  } 
}


int fillpupil(cuFloatComplex *d_odata,float *mask,int nx, int ny, int Nx)
{

  int nnx = nx + BLOCK_SZ - nx%BLOCK_SZ;
  int nny = ny + BLOCK_SZ - ny%BLOCK_SZ;
  dim3 blocks(nnx/BLOCK_SZ,nny/BLOCK_SZ), threads(BLOCK_SZ,BLOCK_SZ);

  fillpupil_krnl<<<blocks, threads>>>(d_odata, mask, nx,ny,Nx);

  cutilCheckMsg("fillpupil_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

__global__ void abs2_krnl(float *odata, cuFloatComplex *idata, int nx, int ny)
{
  __shared__ cuFloatComplex cache[BLOCK_SZ][BLOCK_SZ];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tid = x + y * nx;

  if ((x < nx) && (y < ny)) {
    cache[threadIdx.x][threadIdx.y] = idata[tid];
  } 

  __syncthreads();

  if ((x < nx) && (y < ny)) {
    odata[tid] = (cache[threadIdx.x][threadIdx.y]).x * (cache[threadIdx.x][threadIdx.y]).x +
      (cache[threadIdx.x][threadIdx.y]).y * (cache[threadIdx.x][threadIdx.y]).y;
  } 
}


int abs2(float *d_odata, cuFloatComplex *d_idata, int nx, int ny)
{

  int nnx = nx + BLOCK_SZ - nx%BLOCK_SZ;
  int nny = ny + BLOCK_SZ - ny%BLOCK_SZ;
  dim3 blocks(nnx/BLOCK_SZ,nny/BLOCK_SZ), threads(BLOCK_SZ,BLOCK_SZ);

  abs2_krnl<<<blocks, threads>>>(d_odata, d_idata, nx, ny);

  cutilCheckMsg("abs2_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

