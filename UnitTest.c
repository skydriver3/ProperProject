#ifndef UT 
#define UT 1

#include "SparseTensor.h"
#include "MultiGrid.h"


int DownScale_UpScale_UnitTest(){
  Vector unit_fine; InstantiateVector(&unit_fine, 9*9); 
  ConstantVector(&unit_fine, 1.); 

  Vector unit_coarse; InstantiateVector(&unit_coarse, 5*5); 
  ConstantVector(&unit_coarse, 1.); 

  Vector table_fine; InstantiateVector(&table_fine, 9*9); 
  for(int j = 0 ; j < 9; j++){
    for(int i = 0; i < 9; i++)
    {
      table_fine.vec[i + (j*9)] = i + (j*9); 
    }
  }

  Vector table_broad; InstantiateVector(&table_broad, 5*5); 
  for(int j = 0 ; j < 5; j++){
    for(int i = 0 ; i < 5 ; i++){
      table_broad.vec[i + (j*5)] = i + (j*5); 
    }
  }
  Vector down_out; InstantiateVector(&down_out, 0); 
  DownScaling(unit_fine,&down_out,table_fine, table_broad, 9); 

  Vector up_out; InstantiateVector(&up_out, 0); 
  UpScaling(unit_coarse, &up_out, table_broad, table_fine, 9); 

  return 0; 
}

#endif