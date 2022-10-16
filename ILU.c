#include "SparseTensor.h"


void ILU(MatrixSparse A_original, Vector b, Vector *out){
    MatrixSparse U; InstantiateMatrixSparse(&U, 0); 
    MatrixSparse L; InstantiateMatrixSparse(&L, 0); 
    Identity(&L, b.size); 
    MatrixSparse A; InstantiateMatrixSparse(&A, 0); 
    MatrixCopy(A_original, &A); 
    for (int k = 0; k < b.size; k++){
        int row_start = A.ia[k]; 
        int row_end = A.ia[k+1]; 
        for(int j = row_start ; j < row_end ; j++ ){
            if(A.ja[j] >= k) 
                MatrixInsert(&U, k, A.ja[j], A.a[j]); 

        }
        double akk = MatrixValue(A, k, k); 
        for(int i = k+1 ; i < b.size ; i++){
            double aik = MatrixValue(A, i, k); 
            if(aik != 0)
                MatrixInsert( &L ,i, k, aik/akk); 
            for (int j = k+1 ; j < b.size; j++){
                double aij = MatrixValue(A, i, j); 
                if( aij != 0) {
                    MatrixInsert(&A, i, j, aij - (aik/akk) * MatrixValue(U, k, j)); 
                }
            }
        }
    }

    LUSolve(L, U, b, out); 
}

void LUSolve(MatrixSparse L, MatrixSparse U, Vector b, Vector *out){
    
    
    //Calculating Ux = (1/L)b
    Vector y; InstantiateVector(&y, b.size); 

    for(int i = 0; i < b.size; i++) {
        int row_start = (L.ia)[i] ; 
        int row_end = (L.ia)[i+1]; 
        (y.vec)[i] = (b.vec)[i]; 
        for( int c = row_start; c < row_end; c++){
            if((L.ja)[c] != i)
                (y.vec)[i] -= L.a[c] * (y.vec)[(L.ja)[c]]; 
        }
        (y.vec)[i] /= MatrixValue(L, i+1, i+1); 
    }

    //Calculating x = (1/U)(1/L)b
    for(int i = y.size - 1; i == 0; i--) {
        int row_start = (U.ia)[i] ; 
        int row_end = (U.ia)[i+1]; 
        (out->vec)[i] = (y.vec)[i]; 
        for( int c = row_start; c < row_end; c++){
            if((U.ja)[c] != i)
                (out->vec)[i] -= U.a[c] * (out->vec)[(U.ja)[c]]; 
        }
        (out->vec)[i] /= MatrixValue(U, i+1, i+1); 
    }

}