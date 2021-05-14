//
//  asii_3d.h
//
//  Created by Renaldo Tenório de Moura Júnior and Carlos Vital dos Santos Júnior on 25/09/18.
//  Copyright © 2015 Renaldo Tenorio de Moura Junior. All rights reserved.
//

#ifndef asii_3d_h
#define asii_3d_h
/* In this function the original cube is splicing in n minicubes, the integral value is obtain with the use of the
 * Simple Monte Carlo Numerical Integration Method. Each cub is adapted for diverses iterations for besting of integral. */

	void ASII_3D(double GRIDRECT[20], int NPTS, int NFATX, int NFATY, int NFATZ, int RINT, double dbreak, int refine, double (*fcn)(double []), double *tgral, double *sd){
	
    int i;
    int j;
    int n;
    int k;
    int a;
    int b;
    int c;
    int id;
    int NT_MINICUBS;
    int NPTSTOTAL;
    int *NFAT;
	int index;
	int xi;
	int xj;
	int xk;
	unsigned int seed;

    
    double integral;
    double intresult;
	double intr;
    double result;
	double error;
    double SUM;
    double volume;
	double R;

    double *VECINT;
	double ***RAND;
	double CUR_POINT[4];
    double **RES;
    double ***VECINT3D;
    double **MINICUBSIZE;
    double VECMGRIDRECT[20];
    double **MINIGRIDRECT;
    double **ADPCUBSIZE;
    double *CUBTV;
    double **FAT_EAXIL;

	time_t wall_time;
    double cpu_time;
    double WALLtime;
	double T0 = 0.0;
	double T1 = 0.0;

	//Allocating the vector NFAT and INIPTS for to use the BUILDMINICUBS3D and NNCUBINT
	NFAT = F_VEC_ALLOC_INT(10);
	
	//Initiazing the variables that will be use to refine the Grid
	result=0.0;
	intr=0.0;
	NPTSTOTAL=0;
	SUM = 0.0;
	R = 0.0;
	index=0;

	//Defining Number of the NFAT in x, y and z 
	NFAT[1] = NFATX;
	NFAT[2] = NFATY;
	NFAT[3] = NFATZ;

	//Calculing the total numbers of subspaces
	NT_MINICUBS= NFAT[1]*NFAT[2]*NFAT[3];
	
	//Alocating the Matrix 3D
	VECINT3D = F_3D_ALLOC(NFAT[1]+1,NFAT[2]+1,NFAT[3]+1);
	RAND = F_3D_ALLOC( NT_MINICUBS+1 ,NPTS+1, 3);
	
	//Alocating the Matrix 2D
	MINIGRIDRECT = F_MTX_ALLOC(NT_MINICUBS+1,7);
	FAT_EAXIL =  F_MTX_ALLOC(6,NFAT[1]+1);
	MINICUBSIZE = F_MTX_ALLOC(6,NFAT[1]+1);
	ADPCUBSIZE = F_MTX_ALLOC(6, NT_MINICUBS+1);
	RES =  F_MTX_ALLOC(3, 3*RINT+1);
    
	//Allocating the Vectors
	CUBTV = F_VEC_ALLOC(4);
	VECINT =  F_VEC_ALLOC(NT_MINICUBS+1);
	

	
	//Atibuing the value 0.0 for all spaces of the vector of integrals and 0.0 for vector of subspaces size
	for (i=1 ; i <= NT_MINICUBS; i++){
	
		VECINT[i] = 0.0;
	
	}
	
	//Initializing the matrix 3D VECINT3D
	for (a=0; a <= NFAT[1]; a++){
		for (b=0; b <= NFAT[2]; b++){
			for (c=0; c <= NFAT[3]; c++){
				
				VECINT3D[a][b][c] = 0.0;
				
			}
		}
	}
	
	//Initializating the matrix RES
	for(j=0 ; j < 3*RINT+1; j++ ){	
		for (i=0 ; i < 3; i++ ){

			RES[i][j] =0.0;

		}
	}
	
	//Initializating the matrix RAND
	for(i=0 ; i < NT_MINICUBS+1 ; i++ ){	
		for (j=0 ; j < NPTS+1; j++ ){

			RAND[i][j][0] =0.0;
			RAND[i][j][1] =0.0;
			RAND[i][j][2] =0.0;

		}
	}

	//Initializing the matrix CUR_POINT
	for (i=0 ; i <= 3 ; i++) {
		
		CUR_POINT[i]= 0.0;
			
	}
	
	//Initializating the vector cubtv
	for (i=1 ; i < 4; i++ ){

		CUBTV[i] = 0.0;

	}

	//Atribuition of values for cubtv
	for (i=1 ; i < 4; i++ ){


		CUBTV[i] = GRIDRECT[3 + i] - GRIDRECT[i];

	}

	//Testing tv size
	//Initializing the matrix of delts, MINICUBSIZE and FAT_EAXIL
	for (i= 0 ; i < 6 ; i++){
		for (n=0 ; n <= NFAT[1] ; n++){
			
			MINICUBSIZE[i][n] = 0.0;
			FAT_EAXIL[i][n] = 0.0;
		}
	}
	
	//Initializing the matrix of delts ADPCUBSIZE 
	for (i= 0 ; i < 6 ; i++){
		for (n=0 ; n <= NFAT[1] ; n++){
			
			ADPCUBSIZE[i][n] = 0.0;

		}
	}

	//If CUBTV[1] has the lager size, it is the first
	if (CUBTV[1] > CUBTV[2]) {
		if (CUBTV[1] > CUBTV[3]) {
			
			MINICUBSIZE[1][0] = CUBTV[1];
			MINICUBSIZE[2][0] = CUBTV[2];
			MINICUBSIZE[3][0] = CUBTV[3];

		}
		else
		
		if (CUBTV[3] > CUBTV[2]) {
			
			MINICUBSIZE[1][0] = CUBTV[1];
			MINICUBSIZE[2][0] = CUBTV[3];
			MINICUBSIZE[3][0] = CUBTV[2];
			
		}
		
	}
	
	//If CUBTV[2] has the lager size, it is the first	
	if (CUBTV[2] > CUBTV[1]) {
		if (CUBTV[2] > CUBTV[3]) {
			
			MINICUBSIZE[1][0] = CUBTV[2];
			MINICUBSIZE[2][0] = CUBTV[1];
			MINICUBSIZE[3][0] = CUBTV[3];
			
		}
		else
		
		if (CUBTV[3] > CUBTV[1]) {
			
			MINICUBSIZE[1][0] = CUBTV[2];
			MINICUBSIZE[2][0] = CUBTV[3];
			MINICUBSIZE[3][0] = CUBTV[1];
			
		}
		
	}
	
	//If CUBTV[3] has the lager size, it is the first
	if (CUBTV[3] > CUBTV[1]) {
		if (CUBTV[3] > CUBTV[2]) { 
			
			MINICUBSIZE[1][0] = CUBTV[3];
			MINICUBSIZE[2][0] = CUBTV[1];
			MINICUBSIZE[3][0] = CUBTV[2];
		
		}
		
		else
		
		if (CUBTV[2] > CUBTV[1]) {
			
			MINICUBSIZE[1][0] = CUBTV[3];
			MINICUBSIZE[2][0] = CUBTV[2];
			MINICUBSIZE[3][0] = CUBTV[1];
	
		}
	}

	//Calculating the mini-cub coordinates
	ASII_SLICE(NFAT, GRIDRECT, MINIGRIDRECT, MINICUBSIZE, NT_MINICUBS, CUBTV, FAT_EAXIL);

	//Initialing variables for to use in integration process 
	integral =0.0;
	intresult = 0.0;
	error = 0.0;
	id = 1;
	SUM=0.0;
	NPTS=NPTS*0.1;   // Beging with 10% of the total points number
	
	printf("\n\n");
	printf("                        Iterac.             Integral             Stand. dev. \n");
	printf("                     --------------------------------------------------------\n");

	//Main loop for integration with MS method to OP LMO
	for (n=1; n <= RINT ; n++){

		//Generating random numbers
		for(i=0 ; i < NT_MINICUBS+1 ; i++ ){	
			for (j=0 ; j < NPTS+1; j++ ){

				RAND[i][j][0] = rand();
				RAND[i][j][1] = rand();
				RAND[i][j][2] = rand(); 
				
			}
		}

		intr= 0.0;
		result = 0.0;

		// Open a parallel statement
		#pragma omp parallel default(none) \
					private (seed, k, i, j, xk, xi, xj, VECMGRIDRECT, intr, CUR_POINT, volume, T1, T0,\
					R, integral) \
					shared(RAND, n, NPTS, NT_MINICUBS, MINIGRIDRECT, VECINT, NFAT,fcn) 
			{

			#pragma omp for schedule(dynamic,1) 
			//Open a parallel statement
			for(i=1; i<= NT_MINICUBS; i++){
			
				integral = 0.0;

				//Passing values of the matrix lines that have the mini-cub coordinates 
				for (j=1; j < 7 ; j++) {
					
					VECMGRIDRECT[j] = MINIGRIDRECT[i][j];
				
				}
			
				//---------------------------------------------------------------------------------------------------------
				
				//Calculating the delta for each axis			
				volume = (VECMGRIDRECT[4] - VECMGRIDRECT[1])*(VECMGRIDRECT[5] - VECMGRIDRECT[2])*(VECMGRIDRECT[6] - VECMGRIDRECT[3]);
				// Doing the uniform grid for integration
				for (k=1 ; k <= NPTS ; k++) {

					R = 0.0;
							
					//Calculating curently point
					CUR_POINT[1] = VECMGRIDRECT[1] + ( RAND[i][k][0]/RAND_MAX)*(VECMGRIDRECT[4]-VECMGRIDRECT[1]);
					CUR_POINT[2] = VECMGRIDRECT[2] + ( RAND[i][k][1]/RAND_MAX)*(VECMGRIDRECT[5]-VECMGRIDRECT[2]);
					CUR_POINT[3] = VECMGRIDRECT[3] + ( RAND[i][k][2]/RAND_MAX)*(VECMGRIDRECT[6]-VECMGRIDRECT[3]);
					
					//Funct Calculation
					R = (*fcn)(CUR_POINT);

					integral += R;

				} /////// END NPTS

				//----------------------------------------------------------------------------------------------------------
				
				VECINT[i] = integral*(volume/NPTS);		

			} //////////// END MINICUB				

		}// END OF PRAGMA PARALLEL
		
		for (i=1 ; i <= NT_MINICUBS; i++){
		
			SUM += VECINT[i];
		
		}

		if( refine <= (2*RINT)/3){
			
			intresult +=  SUM/((double)n);

			RES[1][n] = intresult/(double)n;
			error += RES[2][n] = pow((SUM/(double)n)-RES[1][n], 2);
			*sd = sqrt( error/((double)n) );
			*tgral=SUM/(double)n;


		} else {

			intresult +=  SUM/((double)id);

			RES[1][id] = intresult/(double)id;
			error += RES[2][id] = pow((SUM/(double)id)-RES[1][id], 2);
			*sd = sqrt( error/((double)id) );
			*tgral=SUM/(double)id;
	
			id++;
			
		}

		printf("                         %3d:        %16.8lf        %13.6lf\n", n, *tgral , *sd);

		//Definig index value of the matrix 3d VECINT
		index=1;

		//Passing integral value for matrix 3d
		for (a=1; a <= NFAT[1]; a++){
			for (b=1; b <= NFAT[2]; b++){
				for (c=1; c <= NFAT[3]; c++){
					
					VECINT3D[a][b][c] = VECINT[index];
					index++;
					
				} 
			}
		}

		//Adapting each minicub
		ADPT_GRIDRECT3D( 3, NFAT, GRIDRECT, MINIGRIDRECT, NT_MINICUBS, VECINT3D, FAT_EAXIL,
						MINICUBSIZE, ADPCUBSIZE, VECINT, n, dbreak);

		
		if( refine == 2*RINT/3){
			
			NPTS=10*NPTS;
			result = 0.0;
			intresult = 0.0;
			error=0.0;
			SUM=0.0;
			
			printf("\n                         *****      REFINING INTEGRATION GRID!       *****\n\n");
			
		}

		refine++;
		
	}//END 
  
	
	//matrix 3d
	F_3D_DEALLOC(NFAT[1] + 1,NFAT[2] + 1,NFAT[3] + 1, VECINT3D);
    F_3D_DEALLOC(NT_MINICUBS+1,NPTS+1, 3, RAND);

	// Deallocating matrix
	F_MTX_DEALLOC(NT_MINICUBS+1, 7 , MINIGRIDRECT);
	F_MTX_DEALLOC(6, NFAT[1] + 1, MINICUBSIZE);
	F_MTX_DEALLOC(6, NT_MINICUBS+1, ADPCUBSIZE);
	F_MTX_DEALLOC(6, NFAT[1] + 1, FAT_EAXIL);
	F_MTX_DEALLOC(3, 3*RINT+1, RES);

	// Deallocating vectors
	F_VEC_DEALLOC(CUBTV);
	F_VEC_DEALLOC_INT(NFAT);
	F_VEC_DEALLOC(VECINT);
	

}


#endif /* asii_3d_h */