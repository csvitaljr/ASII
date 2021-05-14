//ASII NUMERICAL INTEGRATION MAIN TEST
#include "lib/ASII-NumInt.h"

int main(){
    

    printf("\n    ========================================================================================== ");
    printf("\n    ==============            Adaptive Subspaces by Integral Importance         ============== ");
    printf("\n    ========================================================================================== ");
    printf("\n        > ASII: call test for integrate a normalized gaussian function: ");

    //Defining some local variables
    int i;
    int NPROC    = 1;
    int NINDEX   = 1;                  // Number of threads that goes to each core in a parallel execution
    
    double integral;
    double *GRIDRECT;
    double intgral;
    double sd;
    
    long idumasii;                     // For random number initialization in main. Used in VEGAS and MISER (integrals.h)
 
    // Initializing the random number generator
    srand(time(NULL));
    idumasii = -1*rand( );
  
    // Setting the number of threads equal to NPROC
    NPROC=2;
    omp_set_num_threads(NPROC);
    //printf( "\nNumber of threads used in parallel computation: %d\n\n", NPROC);


    GRIDRECT = F_VEC_ALLOC(20);
    // First put zero in all vactor GRIDRECT
    for (i = 0 ; i < 20 ; i++){
        
        GRIDRECT[i] = 0.0;
        
    }

    //Generating 3D integration space
    GRIDRECT[1]=-10.0; //Xi
    GRIDRECT[2]=-10.0; //Yi
    GRIDRECT[3]=-10.0; //Zi
    GRIDRECT[4]=10.0; //Xf
    GRIDRECT[5]=10.0; //Yf
    GRIDRECT[6]=10.0; //Zf

    //----------------------------------------------------------------------------------------------------------
    //---------------------------------------- BEGING INTEGRAL ASII-3D -----------------------------------------
    /*
    ASII_3D(GRIDRECT, NPTS , NFATX, NFATY, NFATZ, MXINT, dbreak, refine, (*fcn)(double []),  *tgral  , *sd)*/
    ASII_3D(GRIDRECT, 1000 ,   8  ,   8  ,   8  ,  30  ,  0.5  ,   1   ,    &ASII_TEST    , &integral, &sd);

    //---------------------------------------- END INTEGRAL ASII-3D --------------------------------------------
    //----------------------------------------------------------------------------------------------------------

    printf("\n\n\n        > Final Result of ASII numerical integration: integral:%5.6lf sd:%5.6lf ", integral, sd);
    printf("\n\n              Integral:%5.6lf SD:%5.6lf ", integral, sd);
    printf("\n\n\n    ========================================================================================== ");
    printf("\n\n");


    F_VEC_DEALLOC(GRIDRECT);
    
    //**********************************************************************************************************
    
    return(0);
}









