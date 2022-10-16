#ifndef GAUSS_SEIDEL
#define GAUSS_SEIDEL 1
#include "SparseTensor.h"
#include "GaussSeidel.h"
#include <math.h> 
#include <stdio.h> 

/*
This function works out the Gauss Seidel approximation for a given 
iteration using the forward substitution 
"Matrix Computations (3rd edition)" - Golub & Van Loan 1996, eqn (10.1.3)
*/
int GaussSeidelStep(MatrixSparse A, Vector B, Vector *X ){

    // ROW ITERATION 
    for( int i = 0 ; i < B.size; i++){


        //----------TENSOR PRODUCT--------- 
        double sum = 0; 
        int row_start = (A.ia)[i]; 
        int row_end = (A.ia)[i+1];
        for(int j = row_start ; j < row_end ; j++ ){
            int col = (A.ja)[j]; 
            if ( col != i) sum += (X->vec)[col] * (A.a)[j];     
            if ( isnan(sum) ) 
            {
                printf("ERROR : nan during the GaussSeidelStep"); 
                return 1;  
            }     
        }
        //---------------------------------


        (X->vec)[i] = ((B.vec)[i] - sum)/MatrixValue(A, i, i);
        if (isnan(X->vec[i])) {
            printf("ERROR : nan during the GaussSeidelStep"); 
            return 1 ; 
        }
    }
    return 0; 
}

void GaussSeidel(MatrixSparse A, Vector B, Vector *X, int iter){
    for(int i = 0 ; i <= iter ; i++){
        GaussSeidelStep(A, B, X); 
    }
}

#endif //GAUSS_SEIDEL