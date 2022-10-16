#ifndef SPARSETENSOR 
#define SPARSETENSOR 1
#include "SparseTensor.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <math.h>



double MatrixValue(MatrixSparse A, int row, int col){
    int row_start = (A.ia)[row]; 
    int row_end = (A.ia)[row + 1]; 
    for( int c = row_start ; c < row_end; c++) {
        if ( (A.ja)[c] == col) return (A.a)[c]; 
    }
    return 0.0; 
}


void SparseMatrixVecMul(MatrixSparse A, Vector v, Vector *out) {


    for(int row = 0; row < A.ia_size - 1 ; row++){

        double sum = 0; 
        int row_start = (A.ia)[row];
        int row_end = (A.ia)[row+1];
        for (int j = row_start ; j < row_end ; j++) {
            sum += (A.a)[j] * (v.vec)[(A.ja)[j]];
        }

        (out->vec)[row] = sum; 
    }

}

void VecAdd(Vector v1, Vector v2, Vector *out, double coef1, double coef2) {
    
    for(int row = 0; row < v1.size ; row++) {
        (out->vec)[row] = (coef1 * (v1.vec)[row]) + (coef2 * (v2.vec)[row]); 
    }

}

double Dot(Vector a, Vector b) {
    if(a.size != b.size) printf("ERROR : dimension problem in Dot\n"); 

    double sum =0 ; 
    for(int i = 0; i < a.size; i++) {
        sum += (a.vec)[i] * (b.vec)[i]; 
        if(isnan(sum))
        {
            printf("ERROR : encountered a nan in Dot\n"); 
            break; 
        }
    }

    return sum; 
}

double NormVec(Vector v){
    return sqrt(Dot(v, v)); 
}

void Residual(MatrixSparse A, Vector b, Vector u, Vector *r){
    Vector Au; InstantiateVector(&Au, b.size); 
    SparseMatrixVecMul(A, u, &Au); 
    VecAdd(b, Au, r, 1, -1); 
    FreeVector(&Au); 
}

void Zero(Vector *v) {
    for (int i = 0 ; i < v->size ; i++ ) {
        (v->vec)[i] = 0; 
    }
}

void InstantiateVector(Vector *v, int size){
    if(size != 0)
        (v->vec) = (double *) malloc(size * sizeof(double)); 
    (v->size) = size; 
}

void InstantiateMatrixSparse(MatrixSparse *A, int dim){
    (A->ia_size) = 0; 
    (A->ja_size) = 0; 
}

/* 
DEPRECATED use MatrixInsert
row and col start at 1*/ 
int AddElementMatrixSparse(int row, int col, double value, MatrixSparse *A) {
    if(A->ia_size == 0){
        (A->ia) = (int *) malloc(sizeof(int)); 
        (A->ja) = (int *) malloc(sizeof(int)); 
        (A->a) = (double *) malloc((row + 1) * sizeof(double)); 
        (A->ia_size) = row + 1 ; 
        (A->ja_size)++; 
        for(int i = 0 ; i < row; i++){
            (A->ia)[i] = 0 ; 
        }
        (A->ia)[row] = 0; 
        (A->ja)[0] = col - 1; 
        (A->a)[0] = value; 
        return 0; 
    }
    else {
        // Can't insert rows, only append or add to the last one (careful length of ia = number of row + 1) AND Can't add to a column anterior to the last one in the last row
        if(row < (A->ia_size) - 1 || (row == ((A->ia_size) - 1) && col - 1 < (A->ja)[(A->ja_size)-1])) return 1; 

        if(row >= A->ia_size) A->ia = (int *) realloc(A->ia, (row + 1) * sizeof(int)); 
        A->ja = (int *) realloc(A->ja, ((A->ja_size) + 1) * sizeof(int)); 
        A->a = (double *) realloc(A->a, ((A->ja_size) + 1) * sizeof(double)); 

        for(int i = A->ia_size - 1; i<row-1 ; i++ ) {
            (A->ia)[i] = (A->ia)[(A->ia_size) - 2]; 
        }
        (A->ia)[row - 1] = (A->ja_size); 
        (A->ia)[row] = (A->ja_size) + 1; 
        (A->ja)[A->ja_size] = col - 1; 
        (A->a)[A->ja_size] = value; 
        (A->ja_size)++; 
        A->ia_size = row; 
        return 0; 
    }
}

void MatrixCopy(MatrixSparse A, MatrixSparse *out) {
    if(out->ia_size != A.ia_size)
        out->ia = (int *)realloc(out->ia, A.ia_size); 
    if(out->ja_size !=  A.ja_size){
        out->ja = (int *)realloc(out->ja, A.ja_size); 
        out->a = (double *)realloc(out->a, A.ja_size); 
    } 
    memcpy(out->ia, A.ia, A.ia_size * sizeof(int)); 
    memcpy(out->ja, A.ja, A.ja_size * sizeof(int)); 
    memcpy(out->a, A.a, A.ja_size * sizeof(double)); 
    out->ia_size = A.ia_size; 
    out->ja_size = A.ja_size; 
}

int VectorCopy(Vector *dest, Vector src){
    if(dest->size != src.size){
        printf("ERROR: copying vectors didn't work due to size mismatch %d vs %d", dest->size, src.size); 
        return 1; 
    }
    
    memcpy(dest->vec, src.vec, src.size * sizeof(double)); 
    return 0; 
}

/* 
    DEPRECATED Use MatrixInsert
    Can only assign to elements that already exist, thus changing already existing values
   zeros can only stay zeros.

   row and col starts at 1 
   
   return 0 if succesful, 1 otherwise*/
int MatrixAssign(MatrixSparse *A, int row, int col, double value){
    int row_start = (A->ia)[row - 1]; 
    int row_end = (A->ia)[row]; 
    for( int c = row_start ; c < row_end; c++) {
        if ( (A->ja)[c] == col -1 ){
            (A->a)[c] = value; 
            return 0; 
        }
    }
    return 1; 
}

/* 
row and col start at 1 */ 
int MatrixInsert(MatrixSparse *A, int row, int col, double value){
    if(A->ia_size == 0){
        (A->ia) = (int *) malloc((row + 2) *sizeof(int)); 
        (A->ja) = (int *) malloc(sizeof(int)); 
        (A->a) = (double *) malloc(sizeof(double)); 
        (A->ia_size) = row + 2 ; 
        (A->ja_size)++; 
        for(int i = 0 ; i < row + 1; i++){
            (A->ia)[i] = 0 ; 
        }
        (A->ia)[row + 1] = A->ja_size; 
        (A->ja)[0] = col; 
        (A->a)[0] = value; 
        return 0; 
    }
    else if(row + 2 < (A->ia_size) || ( row + 2 == (A->ia_size) && col <= (A->ja)[(A->ja_size)-1])){
        int row_start = (A->ia)[row]; 
        int row_end = (A->ia)[row+1]; 
        for(int j = row_start; j < row_end ; j++) {
            if((A->ja)[j] < col && (j == (row_end-1) || (j != (row_end-1) && (A->ja)[j+1] > col))) {
                //update ja 
                A->ja = (int *) realloc(A->ja, ((A->ja_size) + 1) * sizeof(int) ); 
                A->a = (double *) realloc(A->a, (((A->ja_size)+1) * sizeof(double)));
                for(int k = (A->ja_size); k > j + 1 ; k--) {
                    (A->ja)[k] = (A->ja)[k-1];
                    (A->a)[k] = (A->a)[k-1];
                } 
                (A->ja)[j+1] = col; 
                (A->a)[j+1] = value; 
                //upadate ia
                for(int i = row + 1 ; i < A->ia_size; i++) (A->ia)[i]++; 
            }
            else if (((A->ja)[j] > col && j==row_start)){
                //update ja 
                A->ja = (int *) realloc(A->ja, ((A->ja_size) + 1) * sizeof(int) ); 
                A->a = (double *) realloc(A->a, ((A->ja_size) + 1) * sizeof(double)); 
                A->ja_size++;
                for(int k = (A->ja_size)-1; k > j ; k--) {
                    (A->ja)[k] = (A->ja)[k-1];
                    (A->a)[k] = (A->a)[k-1]; 
                }

                (A->ja)[j] = col; 
                (A->a)[j] = value;
                //upadate ia
                for(int i = row + 1 ; i < A->ia_size; i++) (A->ia)[i]++;                 
            }
            else if((A->ja)[j] == col ) 
                (A->a)[j] = value; 
        }
        if(row_start == row_end){
            for(int k = row + 1 ; k < A->ia_size; k++) {
                (A->ia)[k]++ ;
            } 
            A->ja = (int *) realloc((A->ja),((A->ja_size) + 1) * sizeof(int));
            A->a = (double *) realloc((A->a), ((A->ja_size) + 1) * sizeof(double));  
            A->ja_size++; 
            for(int i = (A->ja_size) - 1; i < (A->ia)[row]; i--){
                A->ja[i] = A->ja[i-1]; 
                A->a[i] = A->a[i-1];
            }
            (A->ja)[(A->ia)[row]] = col; 
            (A->a)[(A->ia)[row]] = value; 
        }
        return 0;
    }
    else{
        A->ia = (int *) realloc(A->ia, (row + 2) * sizeof(int)); 
        A->ja = (int *) realloc(A->ja, ((A->ja_size) + 1) * sizeof(int)); 
        A->a = (double *) realloc(A->a, ((A->ja_size) + 1) * sizeof(double)); 
        A->ja_size++; 

        for(int i = A->ia_size - 1; i<row ; i++ ) {
            (A->ia)[i] = (A->ja_size)-1; 
        }
        (A->ia)[row] = (row + 2) != A->ia_size ? (A->ja_size) - 1 : A->ia[row]; 
        (A->ia)[row+1] = (A->ja_size); 
        (A->ja)[(A->ja_size)-1] = col ; 
        (A->a)[(A->ja_size)-1] = value; 
        A->ia_size = row + 2;
        return 0;          
    }
}

void VectorAppend(Vector *v, double value){
    if(v->size != 0 )
        v->vec = (double *) realloc((v->vec), sizeof(double) * ((v->size) + 1)); 
    else 
        v->vec = (double *) malloc(sizeof(double)); 

    (v->vec)[v->size] = value ; 
    (v->size)++ ; 
}

void SaveVector(const char *filepath, Vector v, Vector table, const int steps){
    FILE *fp; 
    fp = fopen(filepath, "a"); 
    
    if(fp == NULL) {
        if(printf("couldn't open file : %s \n", filepath) <= 0)
        {
            printf("printing did not work \n"); 
        } 
    }
    
    // char *str = (char *) malloc(v.size); 
    for(int i = 0 ; i < table.size; i++){
        fprintf(fp, "%f %d %d \n", v.vec[(int)(table.vec[i])], i % steps, (int)floor(i/steps)) ;  
    }
    fprintf(fp, "\n"); 
    fclose(fp); 
}

void Map2Table(Vector map, Vector *table){
    for(int i = 0 ; i< map.size; i++){
        (table->vec)[(int)((map.vec)[i]) + i] = i ; 
        for(int k = (map.vec)[i] + i + 1; k < (map.vec)[i+1] + i + 1; k++){
            (table->vec)[k] = -1 ; 
        }
    }
}

void ConstantVector(Vector *v, const double val){
    for( int i = 0 ; i < v->size; i++){
        (v->vec)[i] = val ; 
    }
}

int TableInverse(Vector icoor, Vector jcoor, Vector *table, int step){
    int last_ix = -1 ; 
    int ix = 0 ; 
    if( table->size != 0) return 1; 

    for(int k = 0 ; k < icoor.size; k++){
        ix = icoor.vec[k] + jcoor.vec[k] * step ; 

        if (ix - last_ix > 1) {
            for(int i = 0 ; i < ix - last_ix - 1 ; i++){
                VectorAppend(table, -1); 
            }
        }

        VectorAppend(table, k); 
        last_ix = ix ;
    } 

    if( ix < (step - 1) + (step - 1) * step){
        for( int i = 0 ; i < ((step-1) * (step+1)) - last_ix; i++){
            VectorAppend(table, -1); 
        }
    }

    return 0; 
}

void FreeVector(Vector *v){
    free(v->vec);
    v->size = 0 ; 
}

void FreeMatrix(MatrixSparse *A){
    free(A->ja); 
    free(A->a); 
    free(A->ia); 
    A->ia_size = 0 ; 
    A->ja_size = 0; 
}

#endif //SPARSETENSOR