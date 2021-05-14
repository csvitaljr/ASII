//
//  asii_6d.h
//
//  Created by Renaldo Tenorio de Moura Junior and Carlos Vital dos Santos Junior on 25/09/18.
//  Copyright (c) 2015 Renaldo Tenorio de Moura Junior. All rights reserved.
//

#ifndef asii_6d_h
#define asii_6d_h
/* In this function the original cube is splicing in n minicubes, the integral value is obtain with the use of the
 * Simple Monte Carlo Numerical Integration Method. Each cub is adapted for diverses iterations for besting of integral. */

	void ASII_6D(double GRIDRECT[20], int NFATX, int NFATY, int NFATZ, int RINT, double dbreaka, double dbreakb, int refine, int NPTS,double (*fcn)(double [])){
	
	int i;
    int j;
    int k;
    int l;
    int n;
    int a;
    int b;
    int c;
    int d;
	int id;
    int *NFAT;
    int NT_MINICUBS;
    int index;
	int NT_MINICUBS2;
	int xi;
	int xj;
	int xk;
	int adpt;

    double integral;
    double sd;
    double result;
    double SUM;
    double SUM2;
    double SDCOUL;
    double volumeA;
	double volumeB;
	double intresult;
	double intr;
	double error;
	double COUL;
	double da_A, db_A, dc_A;
	double da_B, db_B, dc_B;
	double xe, ye, ze, xe2, ye2, ze2;
    double r_e1_e2; 
	double volsubA;
	double volsubB;
	double R1;
    double R2;
	double CUR_POINT_A[4];
	double CUR_POINT_B[4];
	double VECMGRIDRECT_A[10];
    double VECMGRIDRECT_B[10];
	double XYZ1[4];
    double XYZ2[4];	
    
    double *VECINT_A;
    double *VECINT_B;
    double ***VECINT3D_A;
    double ***VECINT3D_B;
	double ***RAND_A;
	double ***RAND_B;
    double **MINICUBSIZE_A;
    double **MINICUBSIZE_B;
    double **MINIGRIDRECT_A;
    double **MINIGRIDRECT_B;
    double **RES;
    double *CUBTV_A;
    double *CUBTV_B;
    double **FAT_EAXIL_A;
    double **FAT_EAXIL_B;
    double *GRID_CUB_A;
    double *GRID_CUB_B;
    double **INTCOUL_BACKUP;
    double **SDCOUL_BACKUP;
	double TOLZRO = 1e-6;

    
	//Allocating the vector NFAT and INIPTS for to use the BUILDMINICUBS3D and NNCUBINT
	NFAT = F_VEC_ALLOC_INT(10);
	volumeA=0.0;  
	volumeB=0.0;   
	intr = 0.0;
	COUL = 0.0;
	SUM = 0.0;
	SUM2 = 0.0; 
	NT_MINICUBS2 = 0;
	R1 = 0.0;
    R2 = 0.0;
    COUL = 0.0;
	adpt=1;

	//Defining Number of the NFAT in x, y and z 
	NFAT[1] = NFATX;
	NFAT[2] = NFATY;
	NFAT[3] = NFATZ;
	
	//Calculating the total number of integration points
	da_A = db_A = dc_A = 0.0;
	da_B = db_B = dc_B = 0.0;
	volsubA = 0.0;
	volsubB = 0.0;

	//Initializing the matrix CUR_POINT and XYZ1,XYZ2
	for (i=0 ; i <= 3 ; i++) {

		CUR_POINT_A[i]= 0.0;
		CUR_POINT_B[i]= 0.0;

	}

	//Calculating the total numbers of mini-cubes
	NT_MINICUBS= NFAT[1]*NFAT[2]*NFAT[3];
	NT_MINICUBS2 = pow(NT_MINICUBS,2);
	
	//Allocating the Matrix for BUILDMINICUBS_3D
	VECINT3D_A = F_3D_ALLOC(NFAT[1]+1,NFAT[2]+1,NFAT[3]+1);
	VECINT3D_B = F_3D_ALLOC(NFAT[1]+1,NFAT[2]+1,NFAT[3]+1);
	RAND_A = F_3D_ALLOC( NT_MINICUBS2+1 ,NPTS+1, 3);
	RAND_B = F_3D_ALLOC( NT_MINICUBS2+1 ,NPTS+1, 3);

	
	//Allocating the Matrix for BUILDMINICUBS_3D
	MINIGRIDRECT_A = F_MTX_ALLOC(NT_MINICUBS+1,7);
	MINIGRIDRECT_B = F_MTX_ALLOC(NT_MINICUBS+1,7);
	FAT_EAXIL_A =  F_MTX_ALLOC(6,NFAT[1]+1);
	FAT_EAXIL_B =  F_MTX_ALLOC(6,NFAT[1]+1);
	MINICUBSIZE_A = F_MTX_ALLOC(6,NFAT[1]+1);
	MINICUBSIZE_B = F_MTX_ALLOC(6,NFAT[1]+1);
    INTCOUL_BACKUP = F_MTX_ALLOC(NT_MINICUBS+2, NT_MINICUBS+2);
    SDCOUL_BACKUP = F_MTX_ALLOC(NT_MINICUBS+1, NT_MINICUBS+1);
    RES =  F_MTX_ALLOC(3, RINT+1);
    
	//Allocating the Vector for to use the VEGAS function in them
	CUBTV_A = F_VEC_ALLOC(4);
	CUBTV_B = F_VEC_ALLOC(4);
	VECINT_A =  F_VEC_ALLOC(NT_MINICUBS+1);
	VECINT_B =  F_VEC_ALLOC(NT_MINICUBS+1);
    GRID_CUB_A = F_VEC_ALLOC(20);
    GRID_CUB_B = F_VEC_ALLOC(20);

    
	//Initializing the matrix 3D VECINT3D_A and VECINT3D_B
	for (a=0; a <= NFAT[1]; a++){
		for (b=0; b <= NFAT[2]; b++){
			for (c=0; c <= NFAT[3]; c++){
				
				VECINT3D_A[a][b][c] = 0.0;
				VECINT3D_B[a][b][c] = 0.0;
				
			}
		}
	}
	
	//Initializing the vectors VECINT_A and VECINT_B
	for (i=1 ; i <= NT_MINICUBS; i++){
	
		VECINT_A[i] = 0.0;
		VECINT_B[i] = 0.0;
	
	}
	
	//Initializing the vector cubtv
	for (i=1 ; i < 4; i++ ){

		CUBTV_A[i] = 0.0;
		CUBTV_B[i] = 0.0;

	}
	
	//Initializing the vector GRID
	for (i=0 ; i < 20; i++ ){

		GRID_CUB_A[i] = 0.0;
		GRID_CUB_B[i] = 0.0;

	}
	
	//Initializating the vector RES
	for(j=0 ; j < RINT+1; j++ ){	
		for (i=0 ; i < 3; i++ ){

			RES[i][j] =0.0;

		}
	}

	//Attribution of values for cubtv
	for (i=1 ; i < 4; i++ ){

		CUBTV_A[i] = GRIDRECT[6 + i] - GRIDRECT[i];
		CUBTV_B[i] = GRIDRECT[9 + i] - GRIDRECT[3 + i];
		//printf("\n TV    CUBTV_A = %15.6f and CUBTV_B = %15.6f\n", CUBTV_A[i], CUBTV_B[i]);
		
	}
	//Attribution of values for GRID
	for (i=1 ; i < 4; i++){
		
		GRID_CUB_A[i] = GRIDRECT[i];
		GRID_CUB_B[i] = GRIDRECT[i+3];
		GRID_CUB_A[i+3] = GRIDRECT[i+6];
		GRID_CUB_B[i+3] = GRIDRECT[i+9];
	
	}
	
	//Initializing the matrix that keep the coulomb integral and sd for each step
	for (i= 0 ; i < NT_MINICUBS ; i++){
		for (j=0 ; j <= NT_MINICUBS ; j++){
			
			INTCOUL_BACKUP[i][j] = 0.0;
			SDCOUL_BACKUP[i][j] = 0.0;

		}
	}
	
	//------------------------------------------------Testing Tv size------------------------------------------------
	//Initializing the matrix of delts, MINICUBSIZE and FAT_EAXIL
	for (i= 0 ; i < 6 ; i++){
		for (n=0 ; n <= NFAT[1] ; n++){
			
			MINICUBSIZE_A[i][n] = 0.0;
			MINICUBSIZE_B[i][n] = 0.0;
			FAT_EAXIL_A[i][n] = 0.0;
			FAT_EAXIL_B[i][n] = 0.0;
		}
	}

	//------------------------------------------------FOR CUBE "A"------------------------------------------------
	
	//If CUBTV[1] has the lager size, it is the first
	if (CUBTV_A[1] > CUBTV_A[2]) {
		if (CUBTV_A[1] > CUBTV_A[3]) {
			
			MINICUBSIZE_A[1][0] = CUBTV_A[1];
			MINICUBSIZE_A[2][0] = CUBTV_A[2];
			MINICUBSIZE_A[3][0] = CUBTV_A[3];

		}
		else
		
		if (CUBTV_A[3] > CUBTV_A[2]) {
			
			MINICUBSIZE_A[1][0] = CUBTV_A[1];
			MINICUBSIZE_A[2][0] = CUBTV_A[3];
			MINICUBSIZE_A[3][0] = CUBTV_A[2];
			
		}
		
	}
	
	//If CUBTV[2] has the lager size, it is the first	
	if (CUBTV_A[2] > CUBTV_A[1]) {
		if (CUBTV_A[2] > CUBTV_A[3]) {
			
			MINICUBSIZE_A[1][0] = CUBTV_A[2];
			MINICUBSIZE_A[2][0] = CUBTV_A[1];
			MINICUBSIZE_A[3][0] = CUBTV_A[3];
			
		}
		else
		
		if (CUBTV_A[3] > CUBTV_A[1]) {
			
			MINICUBSIZE_A[1][0] = CUBTV_A[2];
			MINICUBSIZE_A[2][0] = CUBTV_A[3];
			MINICUBSIZE_A[3][0] = CUBTV_A[1];
			
		}
		
	}
	
	//If CUBTV[3] has the lager size, it is the first
	if (CUBTV_A[3] > CUBTV_A[1]) {
		if (CUBTV_A[3] > CUBTV_A[2]) { 
			
			MINICUBSIZE_A[1][0] = CUBTV_A[3];
			MINICUBSIZE_A[2][0] = CUBTV_A[1];
			MINICUBSIZE_A[3][0] = CUBTV_A[2];
		
		}
		
		else
		
		if (CUBTV_A[2] > CUBTV_A[1]) {
			
			MINICUBSIZE_A[1][0] = CUBTV_A[3];
			MINICUBSIZE_A[2][0] = CUBTV_A[2];
			MINICUBSIZE_A[3][0] = CUBTV_A[1];
	
		}
	}
    
    //Calculating the subspace coordinates
	ASII_SLICE(NFAT, GRID_CUB_A , MINIGRIDRECT_A, MINICUBSIZE_A, NT_MINICUBS, CUBTV_A, FAT_EAXIL_A);
    
	//------------------------------------------------FOR CUBE "B"------------------------------------------------
	
	//If CUBTV[1] has the lager size, it is the first
	if (CUBTV_B[1] > CUBTV_B[2]) {
		if (CUBTV_B[1] > CUBTV_B[3]) {
			
			MINICUBSIZE_B[1][0] = CUBTV_B[1];
			MINICUBSIZE_B[2][0] = CUBTV_B[2];
			MINICUBSIZE_B[3][0] = CUBTV_B[3];

		}
		else
		
		if (CUBTV_B[3] > CUBTV_B[2]) {
			
			MINICUBSIZE_B[1][0] = CUBTV_B[1];
			MINICUBSIZE_B[2][0] = CUBTV_B[3];
			MINICUBSIZE_B[3][0] = CUBTV_B[2];
			
		}
		
	}
	
	//If CUBTV[2] has the lager size, it is the first	
	if (CUBTV_B[2] > CUBTV_B[1]) {
		if (CUBTV_B[2] > CUBTV_B[3]) {
			
			MINICUBSIZE_B[1][0] = CUBTV_B[2];
			MINICUBSIZE_B[2][0] = CUBTV_B[1];
			MINICUBSIZE_B[3][0] = CUBTV_B[3];
			
		}
		else
		
		if (CUBTV_B[3] > CUBTV_B[1]) {
			
			MINICUBSIZE_B[1][0] = CUBTV_B[2];
			MINICUBSIZE_B[2][0] = CUBTV_B[3];
			MINICUBSIZE_B[3][0] = CUBTV_B[1];
			
		}
		
	}
	
	//If CUBTV[3] has the lager size, it is the first
	if (CUBTV_B[3] > CUBTV_B[1]) {
		if (CUBTV_B[3] > CUBTV_B[2]) { 
			
			MINICUBSIZE_B[1][0] = CUBTV_B[3];
			MINICUBSIZE_B[2][0] = CUBTV_B[1];
			MINICUBSIZE_B[3][0] = CUBTV_B[2];
		
		}
		
		else
		
		if (CUBTV_B[2] > CUBTV_B[1]) {
			
			MINICUBSIZE_B[1][0] = CUBTV_B[3];
			MINICUBSIZE_B[2][0] = CUBTV_B[2];
			MINICUBSIZE_B[3][0] = CUBTV_B[1];
	
		}
	}
	
	volumeA = CUBTV_A[1]*CUBTV_A[2]*CUBTV_A[3];
	volumeB = CUBTV_B[1]*CUBTV_B[2]*CUBTV_B[3];
	
	//Calculating the mini-cub coordinates
	ASII_SLICE(NFAT, GRID_CUB_B, MINIGRIDRECT_B, MINICUBSIZE_B, NT_MINICUBS, CUBTV_B, FAT_EAXIL_B);

	printf("\n\n");
	printf("       Iterac.              Integral             Stand. dev.\n");
	printf("    --------------------------------------------------------\n");

	
    //Initialing integration variables
	integral = 0.0;
	intresult = 0.0;
	error = 0.0;
	id = 1;
	result = 0.0;
	COUL = 0.0;
	NPTS = NPTS*0.1;

	//This loop make n iterations for to adapt the mini-cube coordinates
	for (n=1; n <= RINT; n++){

		//Generating random numbers
		for(i=0 ; i < NT_MINICUBS2+1 ; i++ ){	
			for (j=0 ; j < NPTS+1; j++ ){

				RAND_A[i][j][0] = rand();
				RAND_A[i][j][1] = rand();
				RAND_A[i][j][2] = rand();
				RAND_B[i][j][0] = rand();
				RAND_B[i][j][1] = rand();
				RAND_B[i][j][2] = rand();  

			}
		}
		
		//Initialing local integration variables
		SDCOUL = 0.0;
		intr = 0.0;
		sd = 0.0;
		
		// Open a parallel statement
		#pragma omp parallel default(none) \
					private (r_e1_e2,xe,ye,ze,xe2,ye2,ze2, TOLZRO, XYZ1,XYZ2,R1,R2,integral,CUR_POINT_A, CUR_POINT_B, COUL, \
							n, k, l, i, j, VECMGRIDRECT_A, VECMGRIDRECT_B, volsubA,volsubB,da_B, db_B, dc_B, da_A, db_A, dc_A) \
					shared(RAND_A, RAND_B,    NPTS, NT_MINICUBS, MINIGRIDRECT_A, MINIGRIDRECT_B, INTCOUL_BACKUP,fcn)
		{
			#pragma omp for schedule(dynamic,1) collapse(2)
			//Open a parallel statement, running for all minicubs 
			for(i=1; i<= NT_MINICUBS; i++){
				for(k=1; k <= NT_MINICUBS; k++){
					
					integral = 0.0;
					
					//Passing values of the matrix lines that have the mini-cub coordinates 
					for (l=1; l < 7 ; l++) {
						
						VECMGRIDRECT_A[l] = MINIGRIDRECT_A[i][l];
						VECMGRIDRECT_B[l] = MINIGRIDRECT_B[k][l];
					}
					
					//Calculating the delta for each axis, to A CUBE
					da_A = VECMGRIDRECT_A[4] - VECMGRIDRECT_A[1];
					db_A = VECMGRIDRECT_A[5] - VECMGRIDRECT_A[2];
					dc_A = VECMGRIDRECT_A[6] - VECMGRIDRECT_A[3];
					volsubA = da_A*db_A*dc_A;

					//Calculating the delta for each axis, to B CUBE
					da_B = VECMGRIDRECT_B[4] - VECMGRIDRECT_B[1];
					db_B = VECMGRIDRECT_B[5] - VECMGRIDRECT_B[2];
					dc_B = VECMGRIDRECT_B[6] - VECMGRIDRECT_B[3];
					volsubB = da_B*db_B*dc_B;

					for (j=1 ; j <= NPTS ; j++) {
						
						R1=0.0;
						R2=0.0;
						
						CUR_POINT_A[1] = VECMGRIDRECT_A[1] + ( RAND_A[i*k][j][0]/RAND_MAX)*(VECMGRIDRECT_A[4]-VECMGRIDRECT_A[1]);
						CUR_POINT_A[2] = VECMGRIDRECT_A[2] + ( RAND_A[i*k][j][1]/RAND_MAX)*(VECMGRIDRECT_A[5]-VECMGRIDRECT_A[2]);
						CUR_POINT_A[3] = VECMGRIDRECT_A[3] + ( RAND_A[i*k][j][2]/RAND_MAX)*(VECMGRIDRECT_A[6]-VECMGRIDRECT_A[3]);

						CUR_POINT_B[1] = VECMGRIDRECT_B[1] + ( RAND_B[i*k][j][0]/RAND_MAX)*(VECMGRIDRECT_B[4]-VECMGRIDRECT_B[1]);
						CUR_POINT_B[2] = VECMGRIDRECT_B[2] + ( RAND_B[i*k][j][1]/RAND_MAX)*(VECMGRIDRECT_B[5]-VECMGRIDRECT_B[2]);
						CUR_POINT_B[3] = VECMGRIDRECT_B[3] + ( RAND_B[i*k][j][2]/RAND_MAX)*(VECMGRIDRECT_B[6]-VECMGRIDRECT_B[3]);
						
						//Calculating repulsion between a point in subspace_A with all points of the subspace_B
						//----------------------------------------------------------------------------------------------------------
						
						// Passing XYZ coordinates
						XYZ1[0] = 0.0;
						xe = XYZ1[1] = CUR_POINT_A[1];
						ye = XYZ1[2] = CUR_POINT_A[2];
						ze = XYZ1[3] = CUR_POINT_A[3];
						
						XYZ2[0] = 0.0;
						xe2 = XYZ2[1] = CUR_POINT_B[1];
						ye2 = XYZ2[2] = CUR_POINT_B[2];
						ze2 = XYZ2[3] = CUR_POINT_B[3];
						
						
						// Calculating the distance between the two points
						r_e1_e2 = sqrt( (xe - xe2)*(xe - xe2) +
										(ye - ye2)*(ye - ye2) +
										(ze - ze2)*(ze - ze2) );

						if (r_e1_e2 <= TOLZRO) {
							
							r_e1_e2 = TOLZRO;

						}

						//CALL FUNCTION1
						R1 = (*fcn)(CUR_POINT_A);
						
						//CALL FUNCTION2
						R2 = (*fcn)(CUR_POINT_B);

						COUL = ( R1 * R2 ) / r_e1_e2;
													
						//----------------------------------------------------------------------------------------------------------
							
						integral += COUL; 

					}
					//----------------------------------------------------------------------------------------------------------
					
					INTCOUL_BACKUP[i][k] = (  ( (volsubA*volsubB)/((double)NPTS ) ) * integral  );

				}
			}

		}// END OF PRAGMA PARALLEL

		for (i= 1 ; i <= NT_MINICUBS ; i++){
			for (j=1 ; j <= NT_MINICUBS ; j++){
				
				SUM += INTCOUL_BACKUP[i][j];
				
			}
		}
			
		if( refine <= (2*RINT)/3){
			
			intresult +=  SUM/((double)n);

			RES[1][n] = intresult/(double)n;
			error += RES[2][n] = pow(SUM/((double)n)-RES[1][n], 2);
			result = sqrt( error/((double)n) );

			printf("        %3d:        %16.8lf        %13.6lf\n", n, SUM/(double)n , result);


		} else {

			intresult +=  SUM/((double)id);

			RES[1][id] = intresult/(double)id;
			error += RES[2][id] = pow(SUM/((double)id)-RES[1][id], 2);
			result = sqrt( error/((double)id) );

			printf("        %3d:        %16.8lf        %13.6lf\n", n, SUM/(double)id , result);

			id++;
			
		}

		//Definig index value of the matrix 3d VECINT_A and VECINT_B
		index = 1;

		//Passing integral value for matrix 3d to A CUB
		for (a=1; a <= NFAT[1]; a++){
			for (b=1; b <= NFAT[2]; b++){
				for (c=1; c <= NFAT[3]; c++){
				
					for(d=1; d <= NT_MINICUBS; d++){
					
						VECINT_A[index] += INTCOUL_BACKUP[index][d] ;
						VECINT_B[index] += INTCOUL_BACKUP[index][d] ;

					}
					
					VECINT3D_A[a][b][c] = VECINT_A[index];
					VECINT3D_B[a][b][c] = VECINT_B[index];
					
					index++;
					
				} 
			}
		}
		
		//Decision structure for adaptation. If adpt == 1, else adaptation don't will be activated
		if(adpt == 1){
			//Make the adaptation in cube A, only
			ASII_ADPT6D(NFAT, GRID_CUB_A, MINIGRIDRECT_A, NT_MINICUBS, VECINT3D_A, FAT_EAXIL_A,
									MINICUBSIZE_A, VECINT_A, n, INTCOUL_BACKUP, dbreaka);

			//Make the adaptation in cube B, only
			ASII_ADPT6D(NFAT, GRID_CUB_B, MINIGRIDRECT_B, NT_MINICUBS, VECINT3D_B, FAT_EAXIL_B,
									MINICUBSIZE_B, VECINT_B, n, INTCOUL_BACKUP, dbreakb);
		}

		refine++;
		
		if( refine == (2*RINT)/3 + 1){
		
			NPTS=10*NPTS;
			SUM=0.0;
			result=0.0;
			error=0.03;
			intresult = 0.0;
			adpt=0;
			
			//printf("\n        *****      REFINING INTEGRATION GRID!       *****\n\n");
		
		}
		
	}
	
	// Deallocating Matrix 3d
	F_3D_DEALLOC(NFAT[1] + 1,NFAT[2] + 1,NFAT[3] + 1, VECINT3D_A);
	F_3D_DEALLOC(NFAT[1] + 1,NFAT[2] + 1,NFAT[3] + 1, VECINT3D_B);
	F_3D_DEALLOC(NT_MINICUBS2+1,NPTS+1, 3, RAND_A);
	F_3D_DEALLOC(NT_MINICUBS2+1,NPTS+1, 3, RAND_B);

	// Deallocating matrixs
	F_MTX_DEALLOC(NT_MINICUBS+1, 7 , MINIGRIDRECT_A);
	F_MTX_DEALLOC(NT_MINICUBS+1, 7 , MINIGRIDRECT_B);
	F_MTX_DEALLOC(6, NFAT[1] + 1, MINICUBSIZE_A);
	F_MTX_DEALLOC(6, NFAT[1] + 1, MINICUBSIZE_B);
	F_MTX_DEALLOC(6, NFAT[1] + 1, FAT_EAXIL_A);
	F_MTX_DEALLOC(6, NFAT[1] + 1, FAT_EAXIL_B);
	F_MTX_DEALLOC(NT_MINICUBS+2, NT_MINICUBS+2, INTCOUL_BACKUP);
	F_MTX_DEALLOC(NT_MINICUBS+1, NT_MINICUBS+1,  SDCOUL_BACKUP );
	F_MTX_DEALLOC(3, RINT+1, RES);
	

	// Deallocating vectors
	F_VEC_DEALLOC(CUBTV_A);
	F_VEC_DEALLOC(CUBTV_B);
	F_VEC_DEALLOC_INT(NFAT);
	F_VEC_DEALLOC(VECINT_A);
	F_VEC_DEALLOC(VECINT_B);
	F_VEC_DEALLOC(GRID_CUB_A);
	F_VEC_DEALLOC(GRID_CUB_B);
	
}

#endif /* asii_6d_h */