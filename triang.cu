#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"

//Calcula la distancia entre los puntos |x - y|
__device__ double dist(double x0, double y0, double x1, double y1)
{
	return sqrt(pow(x0 - x1, 2.0) + pow(y0 - y1, 2.0));
}

__device__ int Equality(double a, double b, double epsilon)
{
  return fabs(a - b) < epsilon;
}

//verifica si a es mayor a b de acuerdo a un epsion
__device__ int GreaterEqualthan(double a, double b, double epsilon){
	return Equality(a,b,epsilon) || a > b;
}

/* max_edge_index
 * 
 * Retorna el índice k de la arista máxima de un triángulo i, 
 * descrito por los puntos p0p1p2. Será 0 si p0p1 es máxima.
 * Será 1 si p1p2 lo es. Será 2 si p2p0 lo es.
 * */
__device__ int max_edge_index(int i, double *r, int *p){
     double l0;
     double l1;
     double l2;
     
     int p0;
     int p1;
     int p2;
     
     p0 = p[3*i + 0];
     p1 = p[3*i + 1];
     p2 = p[3*i + 2];
     
     l0 = dist(r[2*p0 + 0], r[2*p0 + 1], r[2*p1 + 0], r[2*p1 + 1]);
     l1 = dist(r[2*p1 + 0], r[2*p1 + 1], r[2*p2 + 0], r[2*p2 + 1]);
     l2 = dist(r[2*p2 + 0], r[2*p2 + 1], r[2*p0 + 0], r[2*p0 + 1]);

     double epsion = 0.001f;
 
     if( (GreaterEqualthan(l0,l1,epsion) && GreaterEqualthan(l1,l2,epsion)) || ( GreaterEqualthan(l0,l2,epsion) && GreaterEqualthan(l2,l1,epsion)))
     {
         return 0;
     }
     else if((GreaterEqualthan(l1,l0,epsion) && GreaterEqualthan(l0,l2,epsion)) || ( GreaterEqualthan(l1,l2,epsion) && GreaterEqualthan(l2,l0,epsion)))
     {
         return 1;
     }
     else
     {
         return 2;
     }
} 


/* same_edge
 * 
 * Indica para las aristas {u,v} y {w,x} si son iguales o no.
 * */
 
 __device__ int same_edge(int u, int v, int w, int x)
 {
     return (u == w && v == x) || (u == x && v == w);
 }


/* get_edge_index
 * 
 * Entrega el índice de la arista {u,v} respecto del triángulo i.
 * */

 __device__  int get_edge_index(int u, int v, int i, int *p)
 {
     int p0;
     int p1;
     int p2;
     
     p0 = p[3*i + 0];
     p1 = p[3*i + 1];
     p2 = p[3*i + 2];
     
     if(same_edge(u, v, p0, p1))
     {
         return 0;
     }
     else if(same_edge(u, v, p1, p2))
     {
         return 1;
     }
     else if(same_edge(u, v, p2, p0))
     {
         return 2;
     }/*
     else
     {
         fprintf(stderr, "%s:%d:%s() ** ERROR ** get_edge_index: Arista {%d,%d} no pertenece al triángulo %d.\n", __FILE__,  __LINE__, __func__, u, v, i);
         exit(EXIT_FAILURE);
     }*/
 }

/* is_nomax_nomax
 * 
 * Indica si la arista compartida entre los triángulos i y j
 * es nomáx-nomáx.
 * */

 __device__ int is_nomax_nomax(int i, int j, int *p, int *max)
 {
     int p0i;
     int p1i;
     int p2i;
     int p0j;
     int p1j;
     int p2j;
     
     p0i = p[3*i + 0];
     p1i = p[3*i + 1];
     p2i = p[3*i + 2];
     
     p0j = p[3*j + 0];
     p1j = p[3*j + 1];
     p2j = p[3*j + 2];
     
     int ij;
     int ii;
     
     if(same_edge(p0i, p1i, p0j, p1j))
     {
         ij = get_edge_index(p0j, p1j, j, p);
         ii = 0;
     }
     else if(same_edge(p1i, p2i, p0j, p1j))
     {
         ij = get_edge_index(p0j, p1j, j, p);
         ii = 1;
     }
     else if(same_edge(p2i, p0i, p0j, p1j))
     {
         ij = get_edge_index(p0j, p1j, j, p);
         ii = 2;
     }
     else if(same_edge(p0i, p1i, p1j, p2j))
     {
         ij = get_edge_index(p1j, p2j, j, p);
         ii = 0;
     }
     else if(same_edge(p1i, p2i, p1j, p2j))
     {
         ij = get_edge_index(p1j, p2j, j, p);
         ii = 1;
     }
     else if(same_edge(p2i, p0i, p1j, p2j))
     {
         ij = get_edge_index(p1j, p2j, j, p);
         ii = 2;
     }
     else if(same_edge(p0i, p1i, p2j, p0j))
     {
         ij = get_edge_index(p2j, p0j, j, p);
         ii = 0;
     }
     else if(same_edge(p1i, p2i, p2j, p0j))
     {
         ij = get_edge_index(p2j, p0j, j, p);
         ii = 1;
     }
     else if(same_edge(p2i, p0i, p2j, p0j))
     {
         ij = get_edge_index(p2j, p0j, j, p);
         ii = 2;
     }
     /*
     else
     {
         fprintf(stderr, "** ERROR ** is_nomax_nomax: Problema insperado para triángulos %d y %d.\n", i, j);
         exit(EXIT_FAILURE);
     }*/
     
     return (ij != max[j]) && (ii != max[i]);
 }

 /* is_max_max
 * 
 * Indica si la arista compartida entre los triángulos i y j
 * es máx-máx.
 * */

 __device__ int is_max_max(int i, int j, int *p, int *max)
{
	int p0i;
	int p1i;
	int p2i;
	
	int p0j;
	int p1j;
	int p2j;
	
	p0i = p[3*i + 0];
	p1i = p[3*i + 1];
	p2i = p[3*i + 2];
	
	p0j = p[3*j + 0];
	p1j = p[3*j + 1];
	p2j = p[3*j + 2];
	
	int ij;
	int ii;
	
	if(same_edge(p0i, p1i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 2;
	}
	else
	{
	//	fprintf(stderr, "** ERROR ** is_max_max: Problema insperado para triángulos %d y %d.\n", i, j);
    //exit(EXIT_FAILURE);
	}
	
	return (ij == max[j]) && (ii == max[i]);
}




/* edge_belongs_to
 * 
 * Indica si arista {k,l} pertenece al triángulo i.
 * */

 __device__ int edge_belongs_to(int k, int l, int i, int *p)
 {
     return same_edge(k, l, p[3*i + 0], p[3*i + 1])
                     || same_edge(k, l, p[3*i + 1], p[3*i + 2])
                     || same_edge(k, l, p[3*i + 2], p[3*i + 0]);
 }

__global__ void label_longest_edges(int *cu_max, double *cu_r, int *cu_triangles)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    cu_max[i] = max_edge_index(i,cu_r, cu_triangles);
}

__global__ void label_frontier_edges(int *cu_max, int *cu_triangles, int *cu_adj, int *cu_seed)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    for(int j = 0; j < 3; j++)
    {
        if(cu_adj[3*i +j] < 0 || is_nomax_nomax(i, cu_adj[3*i + j], cu_triangles, cu_max))
            cu_adj[3*i + j] = NO_ADJ;
        

        if(cu_adj[3*i +j] >= 0 && is_max_max(i, cu_adj[3*i + j], cu_triangles, cu_max) && cu_seed[cu_adj[3*i + j]] == FALSE)
            cu_seed[i] = TRUE;
        
    }   
}

