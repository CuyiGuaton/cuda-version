__device__ int count_FrontierEdges(int triangle, int *cu_adj);

__device__ int generate_polygon(int * poly, int * triangles, int * adj, double *r, int i);

__global__ void generate_mesh(int *cu_triangles, int *cu_adj, double *cu_r, int *cu_seed);
