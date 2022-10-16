#ifndef PROB
#define PROB 1
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prob.h"
#include "SparseTensor.h"
 
#define RADIATOR_WIDTH  0.25




// Permet de raccourcir les arrays jusqu'à la longueur voulu
void TrimArray_int(int ** out, int *toTrim, int length){
    *out = (int *)malloc(sizeof(int) * length); 
    for ( int i = 0; i < length ; i++){
        (*out)[i] = toTrim[i]; 
    }
}

// Permet de raccourcir les arrays jusqu'à la longueur voulu
void TrimArray_double(double **out, double *toTrim, int length){
    *out = (double *) malloc(sizeof(double) * length); 
    for( int i = 0 ; i < length; i++){
        (*out)[i] = toTrim[i]; 
    }
}

// Permet de savoir quelles sont les caractéristiques spatial du point (i, j)
enum topology Project9(int i, int j, double h){

    int n = (int)round( 5 / h ); 
//    printf("h = %f n=%d\n",h, n);
    if( n > i && i > 0 && j == 0 )
        return S; 
    else if ( n > j && j > 0 && i == 0)
        return W; 
    else if ( ( (int)round(3/ h) > i && i > 0 && j == n) || (n > i && i > (int)round(3/h) && j == (int)round(2/h)))
        return N; 
    else if ( ((int)round(2/h) > j && j > 0 && i == n)  || (n > j && j > (int)round(2/h) && i == (int)round(3/h)))
        return E;
    else if ( i == 0 && j == 0)
        return S_W_int;
    else if (i == 0 && j == n) 
        return N_W_int; 
    else if ( i == n && j == (int)round(2/h))
        return N_E_int; 
    else if ( i == (int)round(3/h) && j == n ) 
        return N_E_int; 
    else if ( i == n && j == 0)
        return S_E_int; 
    else if ( i > (int)round(3/h) && i <= n && j > (int)round(2/h) && j <= n ) 
        return Out; 
    return Dom; 
}


// Retourne 1 si le point (i, j) fait parti du domaine de la porte 
int Door(int i, int j, double h) {
    if( (int)round(4./h) >= j && j >= (int)round(2.5 / h) && i == 0)
        return 1; 
    else 
        return 0; 
}

// Retourne 1 si le point (i , j) fait parti du domaine de la fenetre 
int Windows(int i, int j, double h){
    // if( j == 0 && i >= (int)round(0.5 / h) && i <= (int)round(2.5 / h))
    //     return 1 ; 
    // else 
    //     return 0; 
    return 0; 
}


/* Fonction rho qui represente les chauffages */
int Rho(int i, int j, double h) {

    //Chauffage pres de la fenetre 
    // if( i >= (int)round(0.5/h) && i <= (int)round(2.5/h) && j >= (int)round((0.5 - (RADIATOR_WIDTH/2))/ h) && j <= (int)round((0.5 + (RADIATOR_WIDTH / 2))/h))
    //     return 200; 

    // Chauffage pres de la porte 
    // if( i >= (int)round(0.5/h) && i <= (int)round(2.5/h) && j >= (int)round((4.5 - (RADIATOR_WIDTH/2))/ h) && j <= (int)round((4.5 + (RADIATOR_WIDTH / 2))/h))
    //     return 100; 

    //pas de chauffage en cette endroit ( i, j )
    return 0; 
}


/* calcule le nombre de point qu'il ne faut pas prendre en compte se trouvant avant le point i, j
*  On parcourt les i avant de parcourir les j 
*/
void OutsidePointMapping(WallDomain walldom,Domain Door, double h, int **map, int n){
    
    int counter = 0 ; 
    
    for( int j = 0 ; j < n ; j++){
        for( int i = 0 ; i < n ; i++){
            (*map)[i + j * n] = counter;  
            // On verifie si le point (i, j) fait partie du domaine des points dont on cherche la temperature
            if(walldom(i, j, h) != Dom || Door(i, j, h) == 1) 
            {
                // s'il n'en fait pas partie 
                counter++ ;
            }     
        }
    }
}

/* 
Walllength = length of a side of the room in meter 
m = number of steps for a side of the room
n = number of elements
nnz = number of nonzero elements
WallDom = function describing the wall boundaries
WinDom = function describing the window boundaries
DoorDom = function describing the door boundaries
Rho = function describing the heat source 
A = matrix of the problem, should be instantiated at 0 (output)
b = rhs vector should be instantiated at 0(output)
icoor = vector containing the i coordinate of each element of the solution vector, should be instantiated at 0 (output)
jcoor = vector containing the j coordinate of each element of the solution vector, should be instantiated at 0 (output)
*/ 
int prob(double WallLength, int m, int * n, int *nnz, WallDomain WallDom, Domain WinDom, Domain DoorDom, Domain Rho, MatrixSparse *A, Vector *b, Vector *icoor, Vector *jcoor  )
{
    double h = WallLength / ((double)(m - 1));
    double invh2 = 1. / (h * h) ; 
    double Tp = 20.;
    // double Tf = 0.;
    double Tw = 0.;
    int currentLine = 0 ; 

    int * OutsidePointMap; 

    InstantiateMatrixSparse(A, 0);
    InstantiateVector(b, 0); 
    OutsidePointMap = (int *)malloc(m * m * sizeof(int)); 


    OutsidePointMapping(WallDom,DoorDom, h, &OutsidePointMap, m); 

    InstantiateVector(icoor, 0); 
    InstantiateVector(jcoor, 0); 


    // on va parcourir tout les points de la "grille" par laquelle on approxime la piece 
    for( int j = 0 ;  j < m; j++){
        for( int i = 0 ; i < m; i++ ) {
            enum topology _walldom = (*WallDom)(i, j, h); // on regarde la caracteristique spatiale de (i, j)

            int _windom = (*WinDom)(i,j,h); // on regarde si il fait partie de la fenetre 
            int ix = i + j * m ; // l'index dans le vecteur des temperatures sans avoir pris en compte les temperatures qui ne font pas parties du domaines

            int inc = 0 ; 

            if ( (*DoorDom)(i, j, h) == 1 || _walldom != Dom){ // Si on se trouve à une porte ou en dehors du domaine
                inc = 1; // on ne rajouteras l'équation de cette temperature dans notre systeme d'équation 
            }
            else { 


                double gamma = 0. ; // le terme ajouté par les conditions aux bords pour les fenetres 
                double beta = 0.; // le terme ajouté au vecteur b par les conditions aux bords pour les fenetres

                // Pour l'équation numéro "currentLine" on prend la taille actuelle de JA
                int row = currentLine ; 

                // Ajout du terme dans l'équation pour la temperature (i, j-1) 
                if((*DoorDom)(i, j-1, h) == 1){
                    beta += Tp * invh2; //Si (i, j-1) fait parti de la porte alors ne rajoute pas de terme dans A mais un terme dans b qui correspond à la temperature de l'autre piece
                }else if( (*WallDom)(i, j-1, h) != Out && (*WallDom)(i, j-1, h) != Dom ){
                    beta += Tw * invh2; //Si le point au dessus est un mur, la temperature est deja connue, on l'ajoute au membre de droite
                }
                else{
                   
                    int col = ix - m - OutsidePointMap[ix - m]  ; // car il ne faut pas prendre en compte les temperatures qui ne font pas partie du domaine    
                    double value = -invh2 * (1.); //Lorsque (i, j+1) fait parti d'un mur le terme pour (i, j-1) est doublé et vice versa    
                    MatrixInsert(A, row, col, value); 
                }
                    
                
                // idem que (i, j-1) mais pour (i-1, j)
                
                if((*DoorDom)(i - 1, j, h) == 1){
                    beta += Tp * invh2; 
                }else if((*WallDom)(i - 1, j, h) != Out && (*WallDom)(i - 1, j, h) != Dom){
                    beta += Tw * invh2; 
                }
                else{
                    int col = ix - 1 - OutsidePointMap[ix - 1]; 
                    double value = -invh2 * (1.); 
                    MatrixInsert(A, row, col, value); 
                }


                // Ajout du terme pour (i, j)
                int col_center = ix - OutsidePointMap[ix] ; 
                double value_center = (4. + gamma) * invh2; 
                MatrixInsert(A, row, col_center, value_center); 

                // idem que pour (i , j-1) mais pour (i+1, j)    
                if((*DoorDom)(i + 1, j, h) == 1){
                    beta += Tp * invh2; 
                }else if ( ((*WallDom)(i + 1, j, h) != Out) && ((*WallDom)(i + 1, j, h) != Dom)){
                    beta += Tw * invh2; 
                }else{
                    int col = ix + 1 - OutsidePointMap[ix + 1];
                    double value = -invh2 * (1.); 
                    MatrixInsert(A, row, col, value); 
                }


                //idem que pour (i, j-1) mais pour (i, j+1)
                if((*DoorDom)(i , j + 1, h) == 1){
                    beta += Tp * invh2; 
                }else if (((*WallDom)(i, j+1, h) != Out) && ((*WallDom)(i, j+1, h) != Dom)){
                    beta += Tw * invh2; 
                }else
                {                                      
                    int col = ix + m - OutsidePointMap[ix + m]; 
                    double value = -invh2 * ( 1.);                   
                    MatrixInsert(A, row, col, value); 
                }


                VectorAppend(b, beta + Rho(i, j, h)); 
                VectorAppend(icoor, i); 
                VectorAppend(jcoor, j); 
            }
            
            if(inc == 0)
                currentLine++; // on passe à la prochaine équation
        }
    }

    *n = currentLine; 
    *nnz = A->ja_size; 
    return 0; 
}

/* Generator of the matrix, rhs and table for the project 9 
table is links every elements of the solution vector to its index value in the grid*/
void Prob9Gen(int step, MatrixSparse *A, Vector *b, Vector *table) {
    int *n = (int *) malloc(sizeof(int)); 
    int *nnz = (int *) malloc(sizeof(int)); 
    Vector icoor; 
    Vector jcoor; 
    prob(5, step, n, nnz, Project9, Windows, Door, Rho, A, b, &icoor, &jcoor); 
    TableInverse(icoor, jcoor, table, step); 
}


#endif //PROB