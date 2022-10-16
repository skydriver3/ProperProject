#include <stdlib.h> 

void ConvertCSR2CSC(int * ia, int * ja, double * a, int ** colptr, int ** row, double ** val, int n){
    
    int currentRowindex =0; 
    int currentColindex = 0 ; 
    for(int column = 0 ; column < n; column++){
        for(int rowj = 0 ; rowj < n ; rowj++){
            for(int k = ia[rowj]; k < ia[rowj + 1]; k++){
                if(ja[k] == column) row[currentRowindex++] = rowj; 
            }
        }
        colptr[currentColindex++] = currentRowindex; 
    }

}
