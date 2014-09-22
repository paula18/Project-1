#include <iostream>
#include <stdio.h>
#include <cuda.h>
#include <device_launch_parameters.h>

#define SIZE 5
#define BLOCK_DIM 5

__global__ void MatrixAddition(float* d_M, float* d_N, float* d_P)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int index = col + row * SIZE;

    if (col < SIZE && row < SIZE)
    {
        d_P[index] = d_M[index] + d_N[index];
    }
}

__global__ void MatrixSubtraction(float* d_M, float* d_N, float* d_P)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int index = col + row * SIZE;

    if (col < SIZE && row < SIZE)
    {
        d_P[index] = d_M[index] - d_N[index];
    }
}
   

__global__ void MatrixMultiplication(float* d_M, float* d_N, float* d_P, int width)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int index = col + row * width;
   
    int value = 0;

    if (col < width && row < width)
    {
        for (int i = 0; i < width; ++i)
        {
            value += d_M[row * width + i] * d_N[i * width + col];
            d_P[index] = value;
        }
   
    }

}

__host__ void MatAddition(float M[SIZE][SIZE], float N[SIZE][SIZE], float P[SIZE][SIZE], int width)
{
    for ( int i = 0; i < width; ++i )
    {
        for ( int j = 0; j < width; ++j )
        {
            P[j][i] = M[j][i] + N[j][i];
        }
    }

}

__host__ void MatSubtraction(float M[SIZE][SIZE], float N[SIZE][SIZE], float P[SIZE][SIZE], int width)
{
    for ( int i = 0; i < width; ++i )
    {
        for ( int j = 0; j < width; ++j )
        {
            P[j][i] = M[j][i] - N[j][i];
        }
    }

}

__host__ void MatMultiplication(float M[SIZE][SIZE], float N[SIZE][SIZE], float P[SIZE][SIZE], int width)
{
    for ( int i = 0; i < width; ++i )
    {
        for ( int j = 0; j < width; ++j )
        {
            for ( int k = 0; k < width; ++k )
            {
                P[j][i] += M[j][k] * N[k][i];
            }
        }
    }
}


int main()
{
    float m[SIZE][SIZE], n[SIZE][SIZE], pa[SIZE][SIZE], ps[SIZE][SIZE], pm[SIZE][SIZE];
    float pha[SIZE][SIZE], phs[SIZE][SIZE], phm[SIZE][SIZE];
    float *d_ma, *d_na, *d_pa;
    float *d_ms, *d_ns, *d_ps;
    float *d_mm, *d_nm, *d_pm;

    int size = SIZE * SIZE * sizeof(float);
   
    for ( int i = 0; i < SIZE; ++i )
    {
        for ( int j = 0; j < SIZE; ++j )
        {
            m[j][i] = j + i * SIZE;
            n[j][i] = j + i * SIZE;
            pa[j][i] = ps[j][i] = pm[j][i] = pha[j][i] = phs[j][i] = phm[j][i] = 0;
        }
    }

    // Memory allocation

    cudaMalloc(( void**) &d_ma, size );
    cudaMalloc(( void**) &d_na, size );
    cudaMalloc(( void**) &d_pa, size );

    cudaMalloc(( void**) &d_ms, size );
    cudaMalloc(( void**) &d_ns, size );
    cudaMalloc(( void**) &d_ps, size );

    cudaMalloc(( void**) &d_mm, size );
    cudaMalloc(( void**) &d_nm, size );
    cudaMalloc(( void**) &d_pm, size );
   
    cudaMemcpy( d_ma, m, size, cudaMemcpyHostToDevice );
    cudaMemcpy( d_na, n, size, cudaMemcpyHostToDevice );

    cudaMemcpy( d_ms, m, size, cudaMemcpyHostToDevice );
    cudaMemcpy( d_ns, n, size, cudaMemcpyHostToDevice );
   
    cudaMemcpy( d_mm, m, size, cudaMemcpyHostToDevice );
    cudaMemcpy( d_nm, n, size, cudaMemcpyHostToDevice );

    dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
    dim3 dimGrid(1, 1);
   
// Device Operations

    MatrixAddition<<<dimGrid, dimBlock>>>(d_ma, d_na, d_pa);

    MatrixSubtraction<<<dimGrid, dimBlock>>>(d_ms, d_ns, d_ps);

    MatrixMultiplication<<<dimGrid, dimBlock>>>(d_mm, d_nm,d_pm, SIZE);
   
    cudaMemcpy( pa, d_pa, size, cudaMemcpyDeviceToHost );

    cudaMemcpy( ps, d_ps, size, cudaMemcpyDeviceToHost );

    cudaMemcpy( pm, d_pm, size, cudaMemcpyDeviceToHost );

// Host Operations
   
    MatAddition(m, n, pha, SIZE);

    MatSubtraction(m, n, phs, SIZE);

    MatMultiplication(m, n, phm, SIZE);
   
   
    for ( int i = 0; i < SIZE; ++i )
    {
        for ( int j = 0; j < SIZE; ++j )
        {
            //std::cout << i << " " << j << " " << pa[j][i] << std::endl;
            //std::cout << i << " " << j << " " << ps[j][i] << std::endl;
            std::cout << i << " " << j << " " << pm[j][i] << std::endl;

            //std::cout << i << " " << j << " " << pha[j][i] << std::endl;
            //std::cout << i << " " << j << " " << phs[j][i] << std::endl;
            //std::cout << i << " " << j << " " << phm[j][i] << std::endl;

        }
    }

    cudaFree(d_ma);
    cudaFree(d_na);
    cudaFree(d_pa);
   
    cudaFree(d_ms);
    cudaFree(d_ns);
    cudaFree(d_ps);

    cudaFree(d_mm);
    cudaFree(d_nm);
    cudaFree(d_pm);

    system("pause");
}