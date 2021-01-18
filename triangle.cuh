__device__ double dist(double x0, double y0, double x1, double y1);
__device__ int GreaterEqualthan(float a, float b, float epsilon);
__device__ int Equality(double a, double b, double epsilon);
__device__ int max_edge_index(int i, double *r, int *p);
__device__ int is_nomax_nomax(int i, int j, int *p, int *max);
__device__ int get_edge_index(int u, int v, int i, int *p);
__device__ int same_edge(int u, int v, int w, int x);
__device__ int is_max_max(int i, int j, int *p, int *max);
__device__ int edge_belongs_to(int k, int l, int i, int *p);


__global__ void label_longest_edges(int *cu_max, double *cu_r, int *cu_triangles, int tnumber);
__global__ void label_frontier_edges(int *cu_max, int *cu_disconnect, int *cu_triangles, int *cu_adj, int tnumber);
__global__ void get_seeds(int *cu_max, int *cu_triangles, int *cu_adj, int *cu_seed, int tnumber);
__global__ void disconnect_edges(int *cu_adj, int* cu_disconnect, int tnumber);
__global__ void initialize_memory(int *cu_seed, int* cu_disconnect, int tnumber);

