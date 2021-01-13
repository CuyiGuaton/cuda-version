#include <cassert>
#include <iostream>

#define N 100

__global__ void inc(int *a) {
    int i = blockIdx.x;
    if (i<N) {
        a[i]= 1;
    }
}

int main() {
    int ha[N], *da;
    cudaMalloc((void **)&da, N*sizeof(int));
    for (int i = 0; i<N; ++i) {
        ha[i] = i;
    }
    cudaMemcpy(da, ha, N*sizeof(int), cudaMemcpyHostToDevice);
    inc<<<N, 1>>>(da);
    cudaMemcpy(ha, da, N*sizeof(int), cudaMemcpyDeviceToHost);
    for (int i = 0; i < N; ++i) {
        
        std::cout<<ha[i]<<" ";
    }
    
    cudaFree(da);
    return 0;
}