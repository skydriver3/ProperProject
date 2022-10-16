#include <math.h> 
#include <stdlib.h> 

//Multiplication d'une matrice sous format CSR et d'un vecteur de taille respective sizexsize et size
double * Mult( int * ja, double *a , int * ia, double *b, int size){
    double *Vec = (double *) malloc(size * sizeof(double)) ; 

    for(int i = 0 ; i < size; i++){ 
        double tmp = 0.; 
        for( int j = ia[i] ; j < ia[i + 1] ; j++ ) { 
            tmp += b[ja[j]] * a[j]; 
        }
        Vec[i] = tmp; 
    }

    return Vec; 
}

// Différence de deux vecteurs de taille size
double * Subs(double * a, double * b, int size) {
    double * vec = (double *) malloc(size * sizeof(double));
    for(int i = 0 ; i < size; i++) {
        vec[i] = a[i] - b[i]; 
    }
    return vec ; 
}

//Calcule la norme d'un vecteur
double Norm(double *v, int size) {
    double norm = 0. ; 
    for(int i = 0 ; i < size; i++){
        norm += v[i] * v[i]; 
    }
    return sqrt(norm); 
}

// Calcule du résidue 
double Residue(int * ia, int * ja, double * a, double * b, double * u, int size){
    return Norm(Subs(Mult(ja, a, ia, u, size), b, size), size) / Norm(b, size); 
}