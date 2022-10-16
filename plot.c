#include <stdio.h>
#include <stdlib.h>


int plot()
{
    //appel de la commande gnuplot avec un script gnuplot 
//    system("gnuplot GNUscript.p"); 
    system("gnuplot -persistent GNUscript.p"); // AN
    return 0; 
}

