#ifndef PROB_H 
#define PROB_H
#include "SparseTensor.h"

// type de fonction qui vont renvoyer 0 si le point ne se trouve pas dans le domaine et 1 sinon
typedef int (* Domain)(int, int, double); 

enum topology {Dom=0, W = 1, E = -1, S = 2, N = -2, S_W_int = 3, N_W_int = -3, N_E_int = 4, S_E_int = -4, Out = 5 };

// type de fonction qui vont renvoyer
// 0 si le point ne fait pas parti du mur,
// 1 s'il fait parti d'un mur faisant face à l'axe des x vers les positifs 
// -1 s'il fait parti d'un mur faisant face à l'axe des x vers les négatifs 
// 2 s'il fait parti d'un mur faisant face à l'axe des y vers les positifs 
// -2 s'il fait parti d'un mur faisant face à l'axe des y vers les négatifs  
// 3 s'il se trouve dans un coin intérieur orienté vers les positifs 
// -3 s'il se trouve dans un coin intérieur orienté vers les négatifs 
typedef enum topology (* WallDomain)(int, int, double) ;

int prob(double WallLength, int m, int * n, int *nnz, WallDomain WallDom, Domain WinDom, Domain DoorDom, Domain Rho, MatrixSparse *A, Vector *b, Vector *icoor, Vector *jcoor  ); 

int Windows(int i, int j, double h); 

int Door(int i, int j, double h) ; 

int Rho(int i, int j, double h); 

enum topology Project9(int i, int j, double h); 

void Prob9Gen(int step, MatrixSparse *A, Vector *b, Vector *table);

#endif //PROB_H