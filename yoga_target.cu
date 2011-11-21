#include <yoga_target.h>

// declare texture reference for 2D float texture
texture<float, 2, cudaReadModeElementType> tex2d;

// declare shared memory arrays
extern __shared__ cuFloatComplex cachec[];
extern __shared__ float cache[];


__global__ void texraytrace_krnl( float* g_odata, int nx,  int ny, float xoff, float yoff) 
{
  // calculate normalized texture coordinates
  unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
  
  if ((x < nx) && (y < ny)) {
    int xref = x + xoff;
    int yref = y + yoff;
    // read from texture and write to global memory
    g_odata[y*nx + x] += tex2D(tex2d,xref+0.5,yref+0.5);
  }
}

int target_texraytrace(float *d_odata,float *d_idata,int nx, int ny,int Nx, int Ny, float xoff, 
		       float yoff, int Ntot, cudaChannelFormatDesc channelDesc, int device)
{
  tex2d.addressMode[0] = cudaAddressModeClamp;
  tex2d.addressMode[1] = cudaAddressModeClamp;
  tex2d.filterMode = cudaFilterModeLinear;
  tex2d.normalized = false;  

  size_t offset;
  cutilSafeCall(cudaBindTexture2D(&offset,tex2d,d_idata,channelDesc,Nx,Ny,sizeof(float)*Nx));

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  
  int nBlocks = deviceProperties.multiProcessorCount*2;

  dim3 dimBlock(nBlocks,nBlocks , 1);
  dim3 dimGrid(nx / dimBlock.x+1, ny / dimBlock.y+1, 1);

  texraytrace_krnl<<< dimGrid, dimBlock, 0 >>>( d_odata, nx,ny ,xoff,yoff);

  cutilCheckMsg("texraytrace_krnl <<<>>> execution failed\n");

  return EXIT_SUCCESS;
}

__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,float xoff, float yoff, 
			      int Nx,int blockSize)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  int tido;

  int iref,jref,tidi;
  int xref = x + xoff;
  int yref = y + xoff;

  float xshift,yshift,wx1,wx2,wy1,wy2;

  if ((x < nx) && (y < ny)) {
    iref = (int)xref;
    jref = (int)yref;
    tidi = iref + jref * Nx;
    cache[threadIdx.x + threadIdx.y * blockSize] =  idata[tidi];
  } 

  __syncthreads();

  if ((x < nx) && (y < ny)) {
    tido = x + y * nx;

    xshift = xref - iref;
    yshift = yref - jref;

    wx1 = (1.0f - xshift);
    wx2 = xshift;
    wy1 = (1.0f - yshift);
    wy2 = yshift;

    odata[tido] += (wx1 * wy1 * cache[threadIdx.x + threadIdx.y * blockSize] +
		    wx2 * wy1 * cache[threadIdx.x + 1 + threadIdx.y * blockSize] +
		    wx1 * wy2 * cache[threadIdx.x + (threadIdx.y+1) * blockSize] +
		    wx2 * wy2 * cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize]);
  } 
}


int target_raytrace(float *d_odata,float *d_idata,int nx, int ny,int Nx, float xoff, float yoff, int device)
{
  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  /*
  int maxThreads = deviceProperties.maxThreadsPerBlock;
  int nBlocks = deviceProperties.multiProcessorCount*2;
  int nThreads = (N + nBlocks -1)/nBlocks;
  
  if (nThreads > maxThreads) {
    nThreads = maxThreads;
    nBlocks = (N + nThreads  -1)/nThreads;
  }
  */
  // FIX ME !!!!!!!!
  int blockSize = 16;

  int nnx = nx + blockSize - nx%blockSize; // find next multiple of BLOCK_SZ
  int nny = ny + blockSize - ny%blockSize;
  dim3 blocks(nnx/blockSize,nny/blockSize), threads(blockSize,blockSize);

  int smemSize = blockSize * blockSize * sizeof(float);
  raytrace_krnl<<<blocks, threads, smemSize>>>(d_odata, d_idata, nx,ny,xoff,yoff,Nx,blockSize);

  cutilCheckMsg("raytrace_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}


__global__ void fillampli_krnl(cuFloatComplex *odata, float *idata, float *mask, int nx, int ny, 
			       int Nx, int blockSize)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tid = x + y * nx;

  if ((x < nx) && (y < ny)) {
    (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  mask[tid]*cosf(idata[tid]);
    (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  mask[tid]*sinf(idata[tid]);
  } else {
    (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  0.0f;
    (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  0.0f;
  }

  __syncthreads();

  tid = x + y * Nx;

  if ((x < nx) && (y < ny)) {
    odata[tid] = cachec[threadIdx.x + threadIdx.y * blockSize];
  } 
}



int fillampli(cuFloatComplex *d_odata,float *d_idata, float *mask,int nx, int ny, int Nx, int device)
{

  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  // FIX ME !!!!!!!!
  int blockSize = 16;

  int nnx = nx + blockSize - nx%blockSize; // find next multiple of BLOCK_SZ
  int nny = ny + blockSize - ny%blockSize;
  dim3 blocks(nnx/blockSize,nny/blockSize), threads(blockSize,blockSize);

  int smemSize = blockSize * blockSize * sizeof(float);

  cutilCheckMsg("fillampli_kernel<<<>>> execution failed\n");
  fillampli_krnl<<<blocks, threads, smemSize>>>(d_odata, d_idata, mask, nx,ny,Nx, blockSize);

  return EXIT_SUCCESS;
}

__global__ void fillpupil_krnl(cuFloatComplex *odata, float *mask, int nx, int ny, int Nx, int blockSize)
{
  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  //int tid = x + y *blockDim.x * gridDim.x;
  int tid = x + y * nx;

  if ((x < nx) && (y < ny)) {
    (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  mask[tid];
    (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  0.0f;
  } else {
    (cachec[threadIdx.x + threadIdx.y * blockSize]).x =  0.0f;
    (cachec[threadIdx.x + threadIdx.y * blockSize]).y =  0.0f;
  }

  __syncthreads();

  tid = x + y * Nx;

  if ((x < nx) && (y < ny)) {
    odata[tid] = cachec[threadIdx.x + threadIdx.y * blockSize];
  } 
}


int fillpupil(cuFloatComplex *d_odata,float *mask,int nx, int ny, int Nx, int device)
{
  struct cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, device);
  // FIX ME !!!!!!!!
  int blockSize = 16;

  int nnx = nx + blockSize - nx%blockSize; // find next multiple of BLOCK_SZ
  int nny = ny + blockSize - ny%blockSize;
  dim3 blocks(nnx/blockSize,nny/blockSize), threads(blockSize,blockSize);

  int smemSize = blockSize * blockSize * sizeof(float);

  fillpupil_krnl<<<blocks, threads,smemSize>>>(d_odata, mask, nx,ny,Nx,blockSize);

  cutilCheckMsg("fillpupil_kernel<<<>>> execution failed\n");
  return EXIT_SUCCESS;
}

/*
__global__ void texraytrace_krnl2( float* g_odata, int nx,  int ny, float *xref, float *yref) 
{
  // calculate normalized texture coordinates
  unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;
  
  if ((x < nx) && (y < ny))
    // read from texture and write to global memory
    g_odata[y*nx + x] += tex2D(tex,xref[x]+0.5,yref[y]+0.5);
}


__global__ void raytrace_krnl(float *odata, float *idata, int nx, int ny,float xoff, float yoff, 
			      int Nx,int blockSize)
{
  extern __shared__ float cache[];

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;
  int tido;

  int iref,jref,tidi;
  int xref = x + xoff;
  int yref = y + xoff;

  float xshift,yshift,wx1,wx2,wy1,wy2;

  if ((x < nx) && (y < ny)) {
    iref = (int)xref;
    jref = (int)yref;
    tidi = iref + jref * Nx;
    cache[threadIdx.x + threadIdx.y * blockSize] =  idata[tidi];
  } 

  if ((x == nx-1) && (y < ny-1)) {
    cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
    if (threadIdx.y == blockSize-1) {
      cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
      cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
    }
  }

  if ((x < nx-1) && (y == ny-1)) {
    cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
    if (threadIdx.x == blockSize-1) {
      cache[threadIdx.x + 1 + threadIdx.y * blockSize]   = idata[tidi+1];
      cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
    }
  }

  if ((x == nx-1) && (y == ny-1)) {
    cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
    cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
    cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
  }

  __syncthreads();
  
  if ((x < nx) && (y < ny)) {
    if ((threadIdx.x == blockSize-1) && (threadIdx.y == blockSize-1)) {
      cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize] = idata[tidi+Nx+1];
      cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
      cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
    } else {  
      if (threadIdx.x == blockSize-1)
	cache[threadIdx.x + 1 + threadIdx.y * blockSize] = idata[tidi+1];
      if (threadIdx.y == blockSize-1)
	cache[threadIdx.x + (threadIdx.y+1) * blockSize] = idata[tidi+Nx];
    }     
  }
  __syncthreads();

  if ((x < nx) && (y < ny)) {
    tido = x + y * nx;

    xshift = xref - iref;
    yshift = yref - jref;

    wx1 = (1.0f - xshift);
    wx2 = xshift;
    wy1 = (1.0f - yshift);
    wy2 = yshift;

    odata[tido] += (wx1 * wy1 * cache[threadIdx.x + threadIdx.y * blockSize] +
		    wx2 * wy1 * cache[threadIdx.x + 1 + threadIdx.y * blockSize] +
		    wx1 * wy2 * cache[threadIdx.x + (threadIdx.y+1) * blockSize] +
		    wx2 * wy2 * cache[threadIdx.x + 1 + (threadIdx.y+1) * blockSize]);
  } 
}



 */
