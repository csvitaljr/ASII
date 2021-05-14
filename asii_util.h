//
//  asii_util.h
//
//  Extracted from mtx_alloc.h of ChemBOS program lib
//  Copyright © 2016 Renaldo Tenorio de Moura Junior. All rights reserved.
//

#ifndef asii_util_h
#define asii_util_h

#define NR_END 1
#define FREE_ARG char*

/*###########################################################################*/
/*##            double **F_MTX_ALLOC(int _lines, int _rows)                ##*/
/*##                                                                       ##*/
/*##             Function that allocate dynamically a matrix.              ##*/
/*###########################################################################*/

double **F_MTX_ALLOC(int _lines, int _cols) {
    
    int i;  // 'for' counter
    
    double **v;  // Pointer to the matrix
    
    // Checking the received parameters
    if(_lines < 1 || _cols < 1) {
        printf ("*** Error: F_MTX_ALLOC Invalid Parameter ***\n");
        return (NULL);
    }
    
    // Allocating the matrix lines
    v =  calloc (_lines, sizeof(double *));	// Vector of pointers '_lines' to double
    
    if(v == NULL) {
        printf("*** Error: F_MTX_ALLOC Insufficient Memory ***\n");
        return(NULL);
    }
    
    // Allocating the matrix rows
    for( i = 0; i < _lines; i++ ) {
        
        v[i] =  calloc (_cols, sizeof(double));	// Vector of pointers '_rows' to double
        if(v[i] == NULL) {
            printf("*** Error: F_MTX_ALLOC Insufficient Memory ***\n");
            return(NULL);
        }
    }
    
    // Returning the matrix pointer
    return (v);
}  // Endo of double **F_MTX_ALLOC(int _lines, int _rows)


/*###########################################################################*/
/*##     double **F_MTX_DEALLOC(int _linhas, int _colunas, double **v)    ##*/
/*##                                                                       ##*/
/*##           Function that desallocate dynamically a matrix.             ##*/
/*###########################################################################*/
double **F_MTX_DEALLOC(int _lines, int _cols, double **v) {
    
    int  i;  // 'for' counter
    
    
    // Checking the received parameters
    if (v == NULL)
        return (NULL);
    
    if (_lines < 1 || _cols < 1) {
        printf ("*** Error: F_MTX_DEALLOC Invalid Parameter ***\n");
        return(v);
    }
    
    // Freeing the lines
    for (i = 0 ; i < _lines ; i++) {
        free (v[i]);
    }
    
    // Freeing the matrix
    free (v);
    
    // Returning a void pointer
    return (NULL);
    
}  // End of double **F_MTX_DESALLOC(int _linhas, int _colunas, double **v)



/*###########################################################################*/
/*##            double **F_3D_ALLOC(int _lines, int _rows, int _deep)      ##*/
/*##                                                                       ##*/
/*##             Function that allocate dynamically a matrix.              ##*/
/*###########################################################################*/

double ***F_3D_ALLOC(int _lines, int _cols, int _deep) {
    
    int i, j;  // 'for' counter
    
    double ***v;  // Pointer to the matrix
    
    // Checking the received parameters
    if (_lines < 1 || _cols < 1 || _deep < 1) {
        printf ("*** Error: F_3D_ALLOC Invalid Parameter ***\n");
        return (NULL);
    }
    
    // Allocating the matrix lines
    
    // Vector of pointers '_lines' to double
    v =  calloc (_lines, sizeof(double *));
    
    if (v == NULL) {
        printf("*** Error: F_3D_ALLOC Insufficient Memory ***\n");
        return(NULL);
    }
    
    // Allocating the matrix cols
    for ( i = 0; i < _lines; i++ ) {
        
        // Vector of pointers '_cols' to double
        v[i] =  calloc (_cols, sizeof(double));
        
        if (v[i] == NULL) {
            printf("*** Error: F_3D_ALLOC Insufficient Memory ***\n");
            return(NULL);
        }
        
        // Allocating the matrix deep
        for ( j = 0; j < _cols; j++ ) {
            
            // Vector of pointers '_deep' to double
            v[i][j] =  calloc (_deep, sizeof(double));
            
            if (v[i][j] == NULL) {
                printf("*** Error: F_3D_ALLOC Insufficient Memory ***\n");
                return(NULL);
            }
            
        }
        
    }
    
    // Returning the matrix pointer
    return (v);
}  // Endo of double **F_3D_ALLOC(int _lines, int _rows, _deep)


/*###############################################################################*/
/*## double **F_3D_DEALLOC(int _linhas, int _colunas, int _deep, double ***v)  ##*/
/*##                                                                           ##*/
/*##           Function that desallocate dynamically a matrix.                 ##*/
/*###############################################################################*/
double ***F_3D_DEALLOC(int _lines, int _cols, int _deep, double ***v) {
    
    int  i, j;  // 'for' counter
    
    
    // Checking the received parameters
    if (v == NULL)
        return (NULL);
    
    
    // Checking the received parameters
    if (_lines < 1 || _cols < 1 || _deep < 1) {
        printf ("*** Error: F_3D_DEALLOC Invalid Parameter ***\n");
        return(v);
    }
    
    // Freeing the lines
    for (i = 0 ; i < _lines ; i++) {
        
        // Freeing rows
        for (j = 0 ; j < _cols ; j++) {
            
            free (v[i][j]);
            
        }
        
        free (v[i]);
    }
    
    // Freeing the matrix
    free (v);
    
    // Returning a void pointer
    return (NULL);
    
}  // End of double **F_MTX_DESALLOC(int _linhas, int _colunas, double ***v)


/*#########################################################################*/
/*#                  int *F_VEC_ALLOC(int _elementos)                     #*/
/*#                                                                       #*/
/*#      Function that allocate dynamically a FLOAT vector with n ints    #*/
/*#########################################################################*/
double *F_VEC_ALLOC(int _elementos)
{
    
    double *v;        // pointer to vector
    
    // Checking the received parameters
    if (_elementos < 1) {
        printf ("*** Error: F_VEC_ALLOC Invalid Parameter => number of elements ***\n");
        return (NULL);
    }
    
    // allocating the vector
    v =  calloc (_elementos, sizeof(double));
    
    if (v == NULL) {
        printf ("*** Error: F_VEC_ALLOC Insufficient Memory ***\n");
        return (NULL);
    }
    
    // Returning the pointer to vector
    return (v);
}


/*#########################################################################*/
/*#                  double *F_VEC_ALLOC(int _elementos)                   #*/
/*#                                                                       #*/
/*#      Function that desallocate dynamically a vector.                  #*/
/*#########################################################################*/
double *F_VEC_DEALLOC( double *v) {
    
    if (v == NULL) return (NULL);
    free(v);        /* Freeing the vector */
    return (NULL);  /* returns NULL */
    
    
}

/*#########################################################################*/
/*#                  int *F_VEC_ALLOC(int _elementos)                     #*/
/*#                                                                       #*/
/*#      Function that allocate dynamically a FLOAT vector with n ints    #*/
/*#########################################################################*/
int *F_VEC_ALLOC_INT(int _elementos)
{
    
    int *v;        // pointer to vector
    
    // Checking the received parameters
    if (_elementos < 1) {
        printf ("*** Error: F_VEC_ALLOC_INT Invalid Parameter => number of elements ***\n");
        return (NULL);
        
    }
    
    
    // allocating the vector
    v =  calloc (_elementos, sizeof(int));
    
    if (v == NULL) {
        printf ("*** Error: F_VEC_ALLOC_INT Insufficient Memory ***\n");
        return (NULL);
    }
    
    // Returning the pointer to vector
    return (v);
    
}


/*#########################################################################*/
/*#                  double *F_VEC_ALLOC(int _elementos)                   #*/
/*#                                                                       #*/
/*#      Function that desallocate dynamically a vector.                  #*/
/*#########################################################################*/
int *F_VEC_DEALLOC_INT( int *v) {
    
    if (v == NULL) return (NULL);
    free(v);        /* Freeing the vector */
    return (NULL);  /* returns NULL */
}

// ###########################################################################################





#endif /* asii_util_h */
