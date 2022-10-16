#ifndef SPARSETENSOR_H
#define SPARSETENSOR_H


typedef struct MatrixSparse {
    int *ia; 
    int *ja; 
    double *a; 
    int ia_size; 
    int ja_size; 
} MatrixSparse; 

typedef struct Vector {
    double *vec; 
    int size; 
} Vector; 


double MatrixValue(MatrixSparse A, int row, int col);

void SparseMatrixVecMul(MatrixSparse A, Vector v, Vector *out); 

void VecAdd(Vector v1, Vector v2, Vector *out, double coef1, double coef2); 

void InstantiateVector(Vector *v, int size);

void FreeVector(Vector *v);

void FreeMatrix(MatrixSparse *A); 

double Dot(Vector a, Vector b); 

double NormVec(Vector v); 

void Residual(MatrixSparse A, Vector b, Vector u, Vector *r) ; 

void Zero(Vector *v);

int AddElementMatrixSparse(int row, int col, double value, MatrixSparse *A);

void InstantiateMatrixSparse(MatrixSparse *A, int dim);

void MatrixCopy(MatrixSparse A, MatrixSparse *out);

int VectorCopy(Vector *dest, Vector src); 

void VectorAppend(Vector *v, double value);

int MatrixAssign(MatrixSparse *A, int row, int col, double value); 

int MatrixInsert(MatrixSparse *A, int row, int col, double value); 

int TableInverse(Vector icoor, Vector jcoor, Vector *table, int step); 

void SaveVector(const char *filepath, Vector v, Vector table, const int steps);

void ConstantVector(Vector *v, const double val); 

#endif //SPARSETENSOR_H