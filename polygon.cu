#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"
#include "triangle.cuh"


__device__ int count_FrontierEdges(int triangle, int *cu_adj){
    int adj_counter = 0;
    int j;
    for(j = 0; j < 3; j++){ 
        if(cu_adj[3*triangle + j] == NO_ADJ){
            adj_counter++;
        }
    }
    return adj_counter;
}

__device__ int count_BarrierEdges(int *poly, int length_poly){
    int count = 0;
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            count++;
    }
    return count;
}

/*Indica si un triangulo contiene al punto endpoint*/
__device__ int is_continuous(int i, int endpoint, int *p ){
	int p0, p1, p2;
	if (i != -1){
		p0 = p[3*i + 0];
		p1 = p[3*i + 1];
		p2 = p[3*i + 2];
				
		if(endpoint == p0){
			return  0; /* indica que está en p0*/
		}else if (endpoint == p1){
			return  1;  /* indica que está en p1*/
		}else if(endpoint == p2){
			return 2;  /* indica que está en p2*/
		}
	}
	return -1;
}

/* get_adjacent_triangle
 * 
 * Retorna el identificador del triángulo que es adyacente al
 * triángulo i, mediante la arista {k,l}.
 * 
 * Si no hay triángulo, retorna NO_ADY (aún si es porque {k,l}
 * es de borde de triangulación).
 * */

 __device__ int get_adjacent_triangle(int i, int k, int l, int *p, int *adj)
 {
     int u;
     int v;
     int w;
     
     /* Comprobar que {k,l} pertenezca al triángulo i. */
     if(!edge_belongs_to(k, l, i, p))
     {
         //fprintf(stderr, "** ERROR ** get_adjacent_triangle: Arista {%d,%d} no pertenece al triángulo %d.\n", k, l, i);
         //exit(EXIT_FAILURE);
         return -1;
     }
     
     u = adj[3*i + 0];
     v = adj[3*i + 1];
     w = adj[3*i + 2];
     
     int index;
     
     if((u != NO_ADJ)  && (same_edge(k, l, p[3*u + 0], p[3*u + 1]) ||
             same_edge(k, l, p[3*u + 1], p[3*u + 2]) || same_edge(k, l, p[3*u + 2], p[3*u + 0])))
     {
         index = u;
     }
     else if((v != NO_ADJ)  && (same_edge(k, l, p[3*v + 0], p[3*v + 1]) ||
                     same_edge(k, l, p[3*v + 1], p[3*v + 2]) || same_edge(k, l, p[3*v + 2], p[3*v + 0])))
     {
         index = v;
     }
     else if((w != NO_ADJ) && (same_edge(k, l, p[3*w + 0], p[3*w + 1]) ||
                     same_edge(k, l, p[3*w + 1], p[3*w + 2]) || same_edge(k, l, p[3*w + 2], p[3*w + 0])))
     {
         index = w;
     }
     else
     {
         /* Ningún triángulo apareció como adyacente a {k,l} desde el i. */
         index = NO_ADJ;
     }
     
     return index;
 }


/* 
	Busca un triangulo adjacente que comparte el mismo endpoint.
	Origen es el triangulo de donde se viene, -1 si se quiere que se pueda devolver a triangulo anterior.
*/
__device__ int get_adjacent_triangle_share_endpoint(int i, int origen, int endpoint, int *p, int *adj){
	int p0 = p[3*i + 0];
	int p1 = p[3*i + 1];
	int p2 = p[3*i + 2];
	
	/* consigue los triangulos adyacentes */
	int i0 = get_adjacent_triangle(i, p0, p1, p, adj);
	int i1 = get_adjacent_triangle(i, p1, p2, p, adj);
	int i2 = get_adjacent_triangle(i, p2, p0, p, adj);

	/*verifica si los triangulos son continuos al endpoint */
	int ic0 = is_continuous(i0 ,endpoint, p);
	int ic1 = is_continuous(i1 ,endpoint, p);
	int ic2 = is_continuous(i2 ,endpoint, p);
	
	//debug_print("FUNCTION i0 ic0 %d %d   || i1 ic1 %d %d || i2 ic2 %d %d  \n", i0, ic0, i1,ic1,  i2,ic2);
	//debug_print("T %d endpoint %d | Triangles %d %d %d | ADJ  %d %d %d\n", i, endpoint, p[3*i + 0], p[3*i + 1], p[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2] );
	if(ic0 != -1 &&  i0 != origen && i0 != -1){ /*Si hay contuinidad y no retrocede al origen */
		return i0;
	}else if(ic1 != -1 && i1 != origen  && i1 != -1){
		return i1;
	}else if(ic2 != -1 &&   i2 != origen  && i2 != -1){
		return i2;
	}
	return -2;
}

__device__ int generate_polygon(int * poly, int * triangles, int * adj, double *r, int i) {
    int ind_poly = 0;
	
	int initial_point = 0;
	int end_point = 0;
	
	int t0;
	int t1;	
	int t2;
    int ind0;
    int ind1;
    int ind2;
	int continuous;
	int k, j, aux;
	int origen;

    int num_FrontierEdges = count_FrontierEdges(i, adj);
    
    /*si tiene 3 se agregan y se corta el ciclo*/
    if (num_FrontierEdges == 3) {
    
        poly[ind_poly] = triangles[3 * i + 0];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 1];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 2];
        ind_poly++;

        //visited[i] = TRUE;
        return ind_poly;
    } else if(num_FrontierEdges == 2) {
    
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            ind0 = 3*i + j;
            ind1 = 3*i + (j+1)%3;
            ind2 = 3*i + (j+2)%3;
            if(adj[ind0] == NO_ADJ && adj[ind1] == NO_ADJ){
                poly[ind_poly] = triangles[ind1];
                ind_poly++;
                poly[ind_poly] = triangles[ind2];
                ind_poly++;

                initial_point = triangles[ind1];
                end_point = triangles[ind0];  
            }
        }
    }else if (num_FrontierEdges == 1){
    
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            if(adj[3*i + j] == NO_ADJ){
                poly[ind_poly] = triangles[3*i + (j+1)%3];
                ind_poly++;
                initial_point = triangles[3*i + (j+1)%3];

                end_point = triangles[3*i + (j+2)%3];  
            }
        }
    }else {
        end_point = triangles[3*i + 0];
        initial_point = triangles[3*i + 0];
    }
    
    
    /*se marca como visitado */
    //visited[i] = TRUE;
    num_FrontierEdges = 0;
    k = i;
    aux = k;
    k = get_adjacent_triangle_share_endpoint(k, k, end_point, triangles, adj); /* cambia el indice */
    continuous = is_continuous(k, end_point, triangles);
    origen = aux;
//        debug_print("k %d origen %d, conti %d\n", k, origen, continuous);
    
    int triangugulo_initial = i;
    while (initial_point != end_point || triangugulo_initial != k) {

        /*se marca el triangulo visto como visitado y se suma al area del poligono */
        
      //  visited[k] = TRUE;
        t0 = adj[3 * k + 0];
        t1 = adj[3 * k + 1];
        t2 = adj[3 * k + 2];

        num_FrontierEdges = count_FrontierEdges(k, adj);
        
        if (num_FrontierEdges == 2 && continuous != -1) {
            /* ///////////////////si tiene 2 frontier edge se agregan a poly //////////////////////////////////// */

            if (t0 == NO_ADJ && t1 == NO_ADJ) {
                /*si endpoint es continua a t0  y t0-t1 son fe*/
                if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];
                }
            } else if (t2 == NO_ADJ && t0 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2-t0 son fe*/
                if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ && t2 == NO_ADJ) {
                /*si endpoint es continua a t1 y t1-t2 son fe*/
                if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                }
            } 

            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, -1, end_point, triangles, adj); /* se le permite volver al triangulo anterior */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;

        } else if (num_FrontierEdges == 1 && continuous != -1) {
            /* ///////////////////si solo se tiene 1 frontier edge //////////////////////////////////// */
            if (t0 == NO_ADJ) {
                /*si endpoint es continua a t0  y t0 es fe*/
                if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                }
            } else if (t2 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2 es fe*/
                if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ) {
                /*si endpoint es continua a t1  y t1 es fe*/
                if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                }
            } 
            /*si es continuo y tiene 1 fe no puede volver, ind si se guarda  o no*/
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        } else {
            /*si no es continuo no puede regresar de donde venía*/
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        }

    }
    
    return ind_poly;
}


__global__ void generate_mesh(int *cu_triangles, int *cu_adj, double *cu_r, int *cu_seed, int *cu_mesh, int tnumber, int *range){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int i_mesh;
    if(i < tnumber && cu_seed[i]){
        int poly[20]; // CAMBIAR POR SHARE MEMORY
        
        int length_poly = generate_polygon(poly, cu_triangles, cu_adj, cu_r, i);
        __syncthreads(); 
        //cu_mesh[i] = length_poly;

        i_mesh = atomicAdd(range, length_poly+1);
        
        ////guardar malla
        cu_mesh[i_mesh] = length_poly;
        i_mesh++;
        for(int k = 0; k <length_poly; k++){
            cu_mesh[i_mesh + k] = poly[k];
        }
    }
    
}