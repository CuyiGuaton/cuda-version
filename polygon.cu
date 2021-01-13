#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"



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
    //debug_print("Generando polinomio con triangulo %d FE %d\n", i, num_FrontierEdges);
    /*si tiene 3 se agregan y se corta el ciclo*/
    if (num_FrontierEdges == 3) {
        //debug_print("T %d Tiene 3 Frontier edge, se guardan así\n", i);
        poly[ind_poly] = triangles[3 * i + 0];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 1];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 2];
        ind_poly++;

        //visited[i] = TRUE;
        return ind_poly;
    } else if(num_FrontierEdges == 2) {
        //debug_print("T %d Tiene 2 Frontier edge, es oreja, se usa como semilla para generar el poly\n", i);
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
        //debug_print("T %d Tiene 1 Frontier edge,se usa como FE initial\n", i);
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
//        //debug_print("k %d origen %d, conti %d\n", k, origen, continuous);
    //debug_print("T_inicial %d | Triangles %d %d %d | ADJ  %d %d %d\n", i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);
    //debug_print("initial_point %d endpoint %d | T_sig %d\n", initial_point, end_point, k);

    int triangugulo_initial = i;
    while (initial_point != end_point || triangugulo_initial != k) {

        /*se marca el triangulo visto como visitado y se suma al area del poligono */
        
      //  visited[k] = TRUE;
        t0 = adj[3 * k + 0];
        t1 = adj[3 * k + 1];
        t2 = adj[3 * k + 2];

        num_FrontierEdges = count_FrontierEdges(k, adj);
        //debug_print("FE %d | origen %d t %d | Triangles %d %d %d | ADJ  %d %d %d\n", num_FrontierEdges, origen, k, triangles[3*k + 0], triangles[3*k + 1], triangles[3*k + 2], adj[3*k + 0], adj[3*k + 1], adj[3*k + 2]);
        if(origen == -2)
            exit(0);
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
            } else {
                //fprintf(stderr, "** ERROR ** Adding 2 fronter edges\n");
                //fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d \n", k, t0, t1, t2, initial_point, end_point);
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
            } else {
                //fprintf(stderr, "** ERROR ** Adding 1 fronter edges\n");
                //fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d conti %d\n", k, t0, t1, t2, initial_point, end_point, continuous);
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


__global__ void generate_mesh(int *cu_triangles, int *cu_adj, double *cu_r, int *cu_seed){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(cu_seed[i] == TRUE){
        int poly[1000];

        int length_poly = generate_polygon(poly, cu_triangles, cu_adj, cu_r, i);
        int num_BE = count_BarrierEdges(poly, length_poly);
    }
    
}