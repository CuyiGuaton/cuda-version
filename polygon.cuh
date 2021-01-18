__device__ int count_FrontierEdges(int triangle, int *cu_adj);
__device__ int get_adjacent_triangle_share_endpoint(int i, int origen, int endpoint, int *p, int *adj);
__device__ int is_continuous(int i, int endpoint, int *p );
__device__ int get_adjacent_triangle(int i, int k, int l, int *p, int *adj);
__device__ int generate_polygon(int * poly, int * triangles, int * adj, double *r, int i);

__global__ void generate_mesh(int *cu_triangles, int *cu_adj, double *cu_r, int *cu_seed, int *cu_mesh, int tnumber, int *range);
