#ifndef CG_H

#define CG_H 1

#include "MultiGrid.h"
#include "SparseTensor.h"

void ConjugateGradient(MatrixGenerator AGen, Vector *out, Preconditioner preconditioner, int iter, int steps, double precond_relaxation); 

#endif // CG_H 