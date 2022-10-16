#ifndef MG_H
#define MG_H 1

#include "SparseTensor.h"

typedef void (*Smoothing)(MatrixSparse, Vector, Vector *, int); 

typedef void (*Restriction)(Vector, Vector *, Vector, Vector, int); 

typedef int (*Prolongation)(Vector, Vector *, Vector, Vector, int); 

typedef void (*MatrixGenerator)(int, MatrixSparse *, Vector *, Vector *) ; 

typedef int (*Solver)(int, int *, int *, double *, double *,  double *); 

typedef void (*Preconditioner)(MatrixGenerator, MatrixSparse ,Vector ,Vector *, Vector , int , double , int); 

void DownScaling(Vector x, Vector *out, Vector table, Vector table_coarse, int step); 
int UpScaling(Vector broad, Vector *fine, Vector table_broad, Vector table_fine, int step_fine);  
void SaveResidual(int level, MatrixSparse A, Vector b, Vector u, Vector *r, Vector table, const int steps); 
void MG(Smoothing presmoothing, Smoothing postsmoothing, Restriction downsampling, Prolongation upscaling, MatrixGenerator AGen, MatrixSparse A, Vector rhs, Solver solver, Vector *u, Vector table, double relaxation, int level, int depth, int initial_step_number, int internal_epochs); 
void MG_Solver(Smoothing presmoothing, Smoothing postsmoothing, Restriction downsampling, Prolongation upscaling, MatrixGenerator AGen, Solver solver, Vector *out, double relaxation, double acc, int epochs, int depth, int initial_step_number, int internal_epochs); 
void ConvergenceFactor(Vector norms); 

void MG_preconditioner(MatrixGenerator AGen, MatrixSparse A,Vector rhs,Vector *out, Vector table, int steps, double relaxation , int epochs); 

#endif //MG_H