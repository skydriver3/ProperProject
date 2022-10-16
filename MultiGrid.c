#ifndef MULTIGRID 
#define MULTIGRID 1
#include "SparseTensor.h"
#include "prob.h"
#include "MultiGrid.h"
#include "GaussSeidel.h"
#include "umfpk.h"
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>

#define PRESMOOTHING_NU 8
#define POSTSMOOTHING_NU 8
#define SAVE_FILE "residuals/MG.dat"
#define SAVE_FILE_CONVERGENCE "residuals/convergence.dat"

#define RESIDUAL 1

/* Shows in the terminal the configuration of the room */
void TableVisualisation(Vector table, int steps){
    printf("\n STEPS : %d \n", steps);
    for(int j = steps - 1 ; j  >= 0 ; j--){
        for(int i = 0 ; i < steps ; i++){
            printf("%c", (table.vec[i + j * steps] == -1? 'X' : '+')); 
        }
        printf("\n"); 
    }
}


/* Full weighting restriction */
void DownScaling(Vector x, Vector *out, Vector table, Vector table_coarse, int step_fine){
    Vector tmp; InstantiateVector(&tmp, 0); 
    int step_broad = ((step_fine - 1)/2.) + 1; 

    for(int j = 0 ; j < step_fine ; j += 2){
        for(int i = 0; i < step_fine; i += 2){
            
            if(table_coarse.vec[(i/2) + (j/2) * step_broad] != -1){
                int center = (table.vec)[i + j*step_fine]; 
                int right = (i+1) < step_fine ?  (table.vec)[(i+1) + j*step_fine] :  -1; 
                int left = (i-1) >= 0 ? (table.vec)[(i-1) + j*step_fine] : -1; 
                int up = (j-1) >= 0 ? (table.vec)[i + (j-1)*step_fine] : -1; 
                int down = (j+1) < step_fine ? (table.vec)[i + (j+1)*step_fine] : -1; 
                int up_right = (j+1) < step_fine && (i+1) < step_fine ? (table.vec)[(i+1) + (j+1)*step_fine] : -1;
                int down_right = (j-1) >=0 && (i+1) < step_fine ? (table.vec)[(i+1) + (j-1)*step_fine] : -1;
                int up_left = (j+1) < step_fine && (i-1) >= 0 ? (table.vec)[(i-1) + (j+1)*step_fine] : -1;
                int down_left = (j-1) >= 0 && (i-1) >= 0 ? (table.vec)[(i-1) + (j-1)*step_fine] : -1;
                double deno = ( (center != -1 ? 1. : 0.) + (left  != -1? 1. : 0.) + (right != -1? 1. : 0.) + (up != -1? 1. : 0.) + (down != -1? 1. : 0.) + (up_right != -1? 1. : 0.) + (up_left != -1 ? 1. : 0.) + (down_left != -1 ? 1. : 0.) + (down_right != -1 ? 1. : 0.)); 
                double val = ((center != -1 ? (x.vec)[center] : 0.)
                            + (left != -1 ? (x.vec)[left] : 0.)
                            + (right != -1? (x.vec)[right] : 0.)
                            + (up != -1 ? (x.vec)[up] : 0.)
                            + (down != -1 ? (x.vec)[down] : 0.)
                            + (up_right != -1 ? (x.vec)[up_right] : 0.)
                            + (down_right != -1 ? (x.vec)[down_right] : 0.)
                            + (up_left != -1 ? (x.vec)[up_left] : 0.)
                            + (down_left != -1 ? (x.vec)[down_left] : 0.))
                            /  deno ; 

                if(isnan(val)){
                    printf("ERROR : nan in DownScaling");
                    break;
                }

                VectorAppend(out, val ); 
            
                VectorAppend(&tmp, out->size - 1); 
            }
            else
            {
                VectorAppend(&tmp, -1); 
            }
                
            
        }
    }
    //TableVisualisation(tmp, ((step_fine - 1)/2.) + 1); 
    FreeVector(&tmp); 
}

/* Prolongation */
int UpScaling(Vector broad, Vector *fine, Vector table_broad, Vector table_fine, int step_fine){
    int step_broad = ((step_fine - 1)/2.) + 1; 
    for(int j = 0 ; j < step_fine ; j++){
        for(int i = 0; i < step_fine ; i++) {
            int centerf = i + j *step_fine; 
            if(table_fine.vec[centerf] != -1){
                if(i%2 == 0 && j%2 == 0 ){
                    int center = table_broad.vec[(int)(i/2.) + (int)(j/2.) * step_broad]; 
                    double val = (center != -1 ? broad.vec[center] : 0.); 
                    if(isnan(val)){
                        printf("ERROR : nan in UpScaling"); 
                        return 1; 
                    }
                    VectorAppend(fine, val); 
                }
                else if (i%2 == 0 && j%2 == 1 )
                {
                    int up = table_broad.vec[(int)(i/2.) + (int)ceil(j/2.) * step_broad]; 
                    int down = table_broad.vec[(int)(i/2.) + (int)floor(j/2.) * step_broad]; 
                    double val = ((up != -1 ? broad.vec[up] : 0.) + (down != -1 ? broad.vec[down] : 0.)); 
                    val /= (val != 0) ? ( (up != -1 ? 1. : 0.) + (down != -1 ? 1. : 0.)) : 1.; 
                    if(isnan(val)){
                        TableVisualisation(table_broad, step_broad); 
                        printf("coor in fine = i : %d | j : %d", i, j); 
                        printf("ERROR : nan in UpScaling"); 
                        return 1; 
                    }
                    VectorAppend(fine, val); 
                }
                else if (i%2 == 1 && j%2 == 0){
                    int right = table_broad.vec[(int)ceil(i/2.) + (int)(j/2.) * step_broad]; 
                    int left = table_broad.vec[(int)floor(i/2.) + (int)(j/2.) * step_broad]; 
                    double val = (right != -1 ? broad.vec[right] : 0.) + (left != -1 ? broad.vec[left] : 0.); 
                    val /= (val != 0) ? ( (left != -1 ? 1. : 0.) + (right != -1 ? 1. : 0.)) : 1.; 
                    if(isnan(val)){
                        TableVisualisation(table_broad, step_broad) ; 
                        printf("coor in fine = i : %d | j : %d", i, j); 
                        printf("ERROR : nan in UpScaling"); 
                        return 1; 
                    }                    
                    VectorAppend(fine, val); 
                }
                else if (i%2 == 1 && j%2 == 1){
                    int up_right = table_broad.vec[(int)ceil(i/2.) + (int)ceil(j/2.) * step_broad]; 
                    int up_left = table_broad.vec[(int)floor(i/2.) + (int)ceil(j/2.) * step_broad]; 
                    int down_left = table_broad.vec[(int)floor(i/2.) + (int)floor(j/2.) * step_broad]; 
                    int down_right = table_broad.vec[(int)ceil(i/2.) + (int)floor(j/2.) * step_broad]; 
                    double val = (up_right != -1 ? broad.vec[up_right] : 0.) + (up_left != -1 ? broad.vec[up_left] : 0.) + (down_left != -1 ? broad.vec[down_left] : 0.) + (down_right != -1 ? broad.vec[down_right] : 0.); 
                    val /= (val != 0) ? ((up_left != -1 ? 1. : 0.) + (up_right != -1 ? 1. : 0.) + (down_right != -1 ? 1. : 0.) + (down_left != -1 ? 1. : 0.) ) : 1.; 
                    if(isnan(val)){
                        printf("ERROR : nan in UpScaling"); 
                        return 1; 
                    }                    
                    VectorAppend(fine, val); 
                }
                
            }
        }
    }
    return 0; 
}



/*Saves the residual*/
void SaveResidual(int level, MatrixSparse A, Vector b, Vector u, Vector *r, Vector table, const int steps){
    if(!level && RESIDUAL){
        Residual(A, b, u, r); 
        SaveVector(SAVE_FILE, *r, table, steps); 
    }
}
/*Saves the norms of the log10 of the residuals at different iterations */
void ConvergenceFactor(Vector norms){
    printf("The convergence factor : %f \n", (log10(norms.vec[norms.size-1]) - log10(norms.vec[norms.size - 2]))); 
    FILE *fp = fopen(SAVE_FILE_CONVERGENCE, "w"); 
    for(int i = 0 ; i < norms.size ; i++){
        fprintf(fp, "%f\n", log10((norms.vec)[i])); 
    }
    fclose(fp) ; 
}

// As described in "Multigrid Methods" - Volker John
void MG(Smoothing presmoothing, Smoothing postsmoothing, Restriction downsampling, Prolongation upscaling, MatrixGenerator AGen, MatrixSparse A, Vector rhs, Solver solver, Vector *u, Vector table, double relaxation, int level, int depth, int initial_step_number, int internal_epochs){
        
    // TableVisualisation(table, initial_step_number); 

    
    //----------------------Presmoothing of the solution------------------------ 
    presmoothing(A, rhs, u, PRESMOOTHING_NU); 

    //----------------------Restriction on the residual-----------------------
    Vector r; InstantiateVector(&r, rhs.size); 
    Residual(A, rhs, *u, &r); 
    if(!level && RESIDUAL){
        SaveVector(SAVE_FILE, r, table, initial_step_number); 
    }
    Vector r_coarse; InstantiateVector(&r_coarse, 0); 
    Vector table_broad; InstantiateVector(&table_broad, 0); 
    MatrixSparse A_coarse; 
    Vector b_coarse; 
    int step_broad = ((initial_step_number - 1)/2.) + 1; 
    AGen(step_broad, &A_coarse, &b_coarse,&table_broad); 
    //TableVisualisation(table_broad, step_broad); 

    downsampling(r, &r_coarse, table, table_broad, initial_step_number); 
    Vector d_coarse; InstantiateVector(&d_coarse, r_coarse.size);  


    //-----------------------------------------------------------------------
    //---------------------------Solve or Go deeper--------------------------
    if(level == depth ){      
        solver(b_coarse.size, A_coarse.ia, A_coarse.ja, A_coarse.a, r_coarse.vec, d_coarse.vec);
    }  
    else {
        Zero(&d_coarse); 
        for(int i = 0; i < internal_epochs ; i++){
            MG(presmoothing, postsmoothing, downsampling, upscaling, AGen, A_coarse, r_coarse, solver, &d_coarse, table_broad, 1,  level + 1, depth, ((initial_step_number - 1)/2.) + 1, internal_epochs); 
        }
    }

    FreeMatrix(&A_coarse); 
    FreeVector(&b_coarse);
    //---------------------------Apply correction----------------------------
    Vector d; InstantiateVector(&d, 0); 
    upscaling(d_coarse, &d, table_broad, table, initial_step_number); 
    
    VecAdd(*u, d, u, 1, relaxation); 

    SaveResidual(level, A, rhs, *u, &r, table, initial_step_number); 
    
    //----------------------------Post-Smoothing-----------------------------
    postsmoothing(A, rhs, u, POSTSMOOTHING_NU); 

    SaveResidual(level, A, rhs, *u, &r, table, initial_step_number); 

    FreeVector(&d); 
    FreeVector(&r_coarse); 
    FreeVector(&r); 
    FreeVector(&table_broad); 
    FreeVector(&d_coarse); 
 

}


void MG_Solver(Smoothing presmoothing, Smoothing postsmoothing, Restriction downsampling, Prolongation upscaling, MatrixGenerator AGen, Solver solver, Vector *out, double relaxation, double acc, int epochs, int depth, int initial_step_number, int internal_epochs){
    Vector previous_residual; InstantiateVector(&previous_residual, 0); 
    Vector residual_norms; InstantiateVector(&residual_norms, 0); 

    //get the problem matrix and rhs
    MatrixSparse A ; 
    Vector b; 
    Vector table; InstantiateVector(&table, 0); 
    AGen(initial_step_number, &A, &b, &table); 

    //instantiate residual
    InstantiateVector(out, b.size); ConstantVector(out, 20.); 
    Vector residual; InstantiateVector(&residual,b.size); 


    for(int j = 0; j < epochs; j++){

        SaveResidual(0, A, b, *out, &residual, table, initial_step_number);//save the residual for the arbitrary approximation
        MG(presmoothing, postsmoothing, downsampling, upscaling, AGen, A, b, solver, out, table, relaxation, 0, depth, initial_step_number, internal_epochs); //run an iteration of the multigrid method

        if(RESIDUAL){
            
            //store the norm of the residual for this iteration
            VectorAppend(&residual_norms, NormVec(residual)); 
            
            // check for the accuracy of the solution 
            if(!previous_residual.size){
                InstantiateVector(&previous_residual, residual.size); 
                Zero(&previous_residual); 
            }
            Vector res_diff; InstantiateVector(&res_diff, residual.size); 
            VecAdd(residual, previous_residual, &res_diff, 1, -1); 
            if(NormVec(res_diff) > acc){
                FreeVector(&res_diff); 
                VectorCopy(&previous_residual, residual); 
            }
            else 
                break; // if enough get out of the loop 
        }
    }


    if(RESIDUAL) 
        ConvergenceFactor(residual_norms); //save the log10 of the norms for the convergence factor
    
    //free pointers
    FreeVector(&residual), FreeVector(&previous_residual); FreeVector(&residual_norms); 
    FreeVector(&table); FreeVector(&b); 
    FreeMatrix(&A); 
}

void MG_preconditioner(MatrixGenerator AGen, MatrixSparse A,Vector rhs,Vector *out, Vector table, int steps, double relaxation , int epochs){
    for(int i = 0 ; i < epochs; i++){
        MG(GaussSeidel, GaussSeidel, DownScaling, UpScaling, AGen, A, rhs, &solve_umfpack, out, table, relaxation, 0, 5, steps, 2); 
    }
}


#endif //MG