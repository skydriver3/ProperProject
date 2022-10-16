#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include "mkl.h"
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "plot.h"
#include "residue.h"
#include "SparseTensor.h"
#include "MultiGrid.h"
#include "GaussSeidel.h"
#include "ConjugateGradient.h"

int main_MG(){

  printf("\n\n ------------------------ SOLUTION --------------------------\n\n"); 

  Vector solution; InstantiateVector(&solution, 0); 
  //------------------- Parameters ---------------------
  double acc = 0.; // accuracy of the solution
  int m = 257; // number of steps along a dimension of the grid
  int epochs = 25; // number of times to execute back to back the iterative solver on the approximation
  int depth = 5; // =0 : 2-grid || >0 : multi-grid
  int internal_epochs = 2; // =1 : V-cycle || =2 : W-cycle || >2 : custom cycle || ==> represent the number of times MG should be used on the coarser grid
  //-------------------- Relaxation --------------------
  double slope = 1.66; 
  double convergence = pow(10, -slope); 
  double lambda_max = 1.01579 ; 
  double lambda_min = lambda_max * ( (1. - convergence)/(1. + convergence) ); 
  double relaxation = 2. / (lambda_max + lambda_min);
  double th_asym_conv_factor = 1. - (relaxation * lambda_min); 
  // relaxation = 1; 
  printf("Relaxation = %f\n", relaxation); 
  printf("Convergence = %f\n", convergence); 
  printf("lambda min = %f\n", lambda_min); 
  printf("the theoretical accelerated asymptotic convergence factor is %f and the corresponding slope is %f \n", th_asym_conv_factor, log10(1/th_asym_conv_factor)); 
  //----------------------- Solver ---------------------
  double t1 = mytimer(); 
  MG_Solver(&GaussSeidel, &GaussSeidel, &DownScaling, &UpScaling, &Prob9Gen, &solve_umfpack, &solution, relaxation,  acc, epochs, depth, m, internal_epochs); 
  double t2 = mytimer();
  printf("\nTime taken to compute solution (CPU): %f sec\n",t2-t1); 
  FreeVector(&solution); 

  printf("\n\n ==================================================================\n\n"); 

  return 0; 
}

int main_CG(){

  printf("\n\n ------------------------ SOLUTION --------------------------\n\n"); 

  Vector solution; InstantiateVector(&solution, 0); 
  //------------------- Parameters ---------------------
  double acc = 0.; // accuracy of the solution
  int m = 257; // number of steps along a dimension of the grid
  int epochs = 30; // number of times to execute back to back the iterative solver on the approximation
  int depth = 5; // =0 : 2-grid || >0 : multi-grid
  int internal_epochs = 2; // =1 : V-cycle || =2 : W-cycle || >2 : custom cycle || ==> represent the number of times MG should be used on the coarser grid
  //-------------------- Relaxation --------------------
  double slope = 1.66; 
  double convergence = pow(10, -slope); 
  double lambda_max = 1.01579 ; 
  double lambda_min = lambda_max * ( (1. - convergence)/(1. + convergence) ); 
  double relaxation = 2. / (lambda_max + lambda_min);
  double th_asym_conv_factor = 1. - (relaxation * lambda_min); 
  double sqrt_K = sqrt(lambda_max / lambda_min); 
  double max_slope = log10(2 * ( (sqrt_K - 1) / (sqrt_K + 1))); 
  relaxation = 1; 
  printf("Relaxation = %f\n", relaxation); 
  printf("Convergence = %f\n", convergence); 
  printf("lambda min = %f\n", lambda_min); 
  printf("the theoretical accelerated asymptotic convergence factor is %f and the corresponding slope is %f \n", th_asym_conv_factor, log10(1/th_asym_conv_factor)); 
  printf("\nThe maximum slope for the logarithm of the residual is %f\n", max_slope);
  //----------------------- Solver ---------------------
  double t1 = mytimer(); 
  ConjugateGradient(&Prob9Gen, &solution, &MG_preconditioner, epochs, m, relaxation); 
  double t2 = mytimer();
  printf("\nTime taken to compute solution (CPU): %f sec\n",t2-t1); 
  FreeVector(&solution); 

  printf("\n\n ==================================================================\n\n"); 

  return 0; 
}

int main(){
  main_CG(); 
}


