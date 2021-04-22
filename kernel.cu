/*Vector add kernel*/
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"

using namespace std;


__global__ void vectorAdd(int* a, int* b, int* c, int n) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if(tid < n)
		c[tid] = a[tid] + b[tid];
}

__global__ void matrixMul(int* m, int* n, int* p, int ns) {
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int column = blockIdx.x * blockDim.x + threadIdx.x;

	int p_sum = 0;
	if(row<ns && column<ns) {
		for (int i = 0; i < ns; i++) {
		p_sum += m[row * ns + i] * n[i * ns + column];
	}
	p[row * ns + column] = p_sum;
}

int main() {
	int n = 1 << 10; 


	//Host pointers
	int* h_a;
	int* h_b;
	int* h_c;

	//Device Pointers
	int* d_a;
	int* d_b;
	int* d_c;

	size_t bytes = n * sizeof(int);

	// Allocate memory on host side
	h_a = (int*)malloc(bytes);
	h_b = (int*)malloc(bytes);
	h_c = (int*)malloc(bytes);

	for (int i = 0; i < n; i++) {
		h_a[i] = rand();
		h_b[i] = rand();
	}

	// Allocate memory on device side
	cudaMalloc(&d_a, bytes);
	cudaMalloc(&d_b, bytes);
	cudaMalloc(&d_c, bytes);

	// Init block and grid size
	int block_size = 1024;
	int grid_size = (int)ceil( (float)n / block_size);
	printf("Grid size if %d\n", grid_size);

	cudaMemcpy(d_a, h_a, bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, bytes, cudaMemcpyHostToDevice);

	vectorAdd<<<grid_size, block_size>>> (d_a, d_b, d_c, n);

	cudaMemcpy(h_c, d_c, bytes, cudaMemcpyDeviceToHost);


	free(h_a);
	free(h_b);
	free(h_c);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	//Host Matrix m, n, p
	int* h_m;
	int* h_n;
	int* h_p;


	//Device Matrix m, n, p
	int* d_m;
	int* d_n;
	int* d_p;

	size_t bytesM = n * n * sizeof(int);

	h_n = (int*)malloc(bytesM);
	h_m = (int*)malloc(bytesM);
	h_p = (int*)malloc(bytesM);

	//Initialize Matrix
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			h_n[i * n + j] = rand() % 1024;
			h_m[i * n + j] = rand() % 1024;
		}
	}

	cudaMalloc(&d_m, bytesM);
	cudaMalloc(&d_n, bytesM);
	cudaMalloc(&d_p, bytesM);

	//Copy data to the device
	cudaMemcpy(d_m, h_m, bytesM, cudaMemcpyHostToDevice);
	cudaMemcpy(d_n, h_n, bytesM, cudaMemcpyHostToDevice);

	int threads_per_block = 16;
	dim3 block_sizeD(threads_per_block, threads_per_block);
	dim3 grid_sizeD(n / block_sizeD.x, n / block_sizeD.y);

	matrixMul << <grid_sizeD, block_sizeD >> > (d_m, d_n, d_p, n);
	cudaMemcpy(h_p, d_p, bytesM, cudaMemcpyDeviceToHost);

	free(h_m);
	free(h_n);
	free(h_p);

	cudaFree(d_m);
	cudaFree(d_n);
	cudaFree(d_p);
	return 0;
}
