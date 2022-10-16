#include "SparseTensor.h"
#include "ConjugateGradient.h"
#include "MultiGrid.h"
#include "prob.h"
#include <stdio.h>
#include <math.h>

#define CG_SAVE_FILE "residuals/CG.dat"
#define CG_CONV_SAVE_FILE "residuals/convergence_CG.dat"

void ConjugateGradient(MatrixGenerator AGen, Vector *out, Preconditioner preconditioner, int iter, int steps, double precond_relaxation) {

    //================ INITIALIZATION ======================

    Vector residual_norms; InstantiateVector(&residual_norms, 0) ; 

    MatrixSparse A; 
    Vector B; 
    Vector table; InstantiateVector(&table, 0);  
    AGen(steps, &A, &B, &table); 

    Vector Real_residual; InstantiateVector(&Real_residual, B.size); 
    
    // -----------Arbitrary first approximation-----------
    InstantiateVector(out, B.size); 
    ConstantVector(out, 20.); 

    Residual(A, B, *out, &Real_residual); 
    VectorAppend(&residual_norms, NormVec(Real_residual)); 
    // -----------Residual from first approximation-------
    Vector r; InstantiateVector(&r, B.size); 
    Residual(A, B, *out, &r); 

    
 
    Vector previous_d; InstantiateVector(&previous_d, r.size);
    double previous_alpha_deno = 0 ; 

    //================== ITERATION ===========================

    for(int i = 0; i < iter ; i++) {
        //----------------Beta------------------------
        Vector precond_r; InstantiateVector(&precond_r, r.size); 
        preconditioner(AGen, A, r, &precond_r, table, steps, precond_relaxation, 5); 
        Vector ABr; InstantiateVector(&ABr, precond_r.size); 
        SparseMatrixVecMul(A, precond_r, &ABr); 
        double beta_num = Dot(previous_d, ABr) ;

        double beta = (i==0)? 0 : -beta_num / previous_alpha_deno; 
        //-----------------Gradient--------------------
        Vector d; InstantiateVector(&d, r.size); 
        VecAdd(precond_r, previous_d, &d, 1, (i==0)?0:beta); 
        if(NormVec(d) == 0){
            break; 
        } 

        //------------------Alpha----------------------
        Vector Ad; InstantiateVector(&Ad, d.size); 
        SparseMatrixVecMul(A, d, &Ad); 
        previous_alpha_deno = Dot(d, Ad); 
        double alpha = Dot(d, r) / previous_alpha_deno; 

        //------------------Update Solution-------------
        VecAdd(*out, d, out, 1, alpha); 

        //------------------Update Residual-------------
        VecAdd(r, Ad, &r, 1, -alpha);

        //------------------Save Residual-----------------
        Residual(A, B, *out, &Real_residual); 
        VectorAppend(&residual_norms, NormVec(Real_residual)); 

        VectorCopy(&previous_d, d); 
        FreeVector(&precond_r); FreeVector(&d); FreeVector(&Ad); 
    }



    //----------------------Save convergence-----------------
    printf("The convergence factor : %f \n", (log10(residual_norms.vec[residual_norms.size-1]) - log10(residual_norms.vec[residual_norms.size - 2]))); 
    FILE *fp = fopen(CG_CONV_SAVE_FILE, "w"); 
    for(int i = 0 ; i < residual_norms.size ; i++){
        fprintf(fp, "%f\n", log10((residual_norms.vec)[i])); 
    }
    fclose(fp) ; 

    //------------------------Freeing Vector-------------------
    FreeVector(&r);
    FreeVector(&previous_d); 
    FreeVector(&residual_norms); 
    FreeVector(&B); 
    FreeMatrix(&A); 
    FreeVector(&Real_residual); 
}

