//
//  asii_adpt6d.h
//
//  Created by Renaldo Tenorio de Moura Junior and Carlos Vital dos Santos Junior on 13/05/21.
//  Copyright © 2021 Renaldo Tenorio de Moura Junior. All rights reserved.
//
//  This function adapt the  coordinates that define the size of each mini-cub for a good parallel balanncig

#ifndef asii_adpt6d_h
#define asii_adpt6d_h

// NFAT is the number of the mini-cub that will be generate
// GRIDRECT3D contains the cub-coordinates

 void ASII_ADPT6D(int NFAT[4], double GRIDRECT3D[20], double **MINIRECTSIZE, int NT_MINICUB, double ***VECINT3D, 
 						double **FAT_EAXIL, double **MINICUBSIZE, double *VECINT, int NIT_INTGR, double **INTCOUL_BACKUP, 
						double dbreak){
	
	//Defining variables here used
	int i;
	int n;
	int j;
	int k;
	int a;
	int b;
	int c;
	int d;
	int l;
	int count;
	int index;
	int freeze;

	double SUMdelt[5];
	double FAT_INTGR[4][NFAT[1]+1];
	double integral;
	double sd;
	double intresult;
	double intsd;
	double VECMGRIDRECT_A[10];
	double VECMGRIDRECT_B[10];
	double SUM;
	double expt;
    double initotal;
	double fatnorm;
	double volume;
	double *VOLUME;
	
	//Calculating the tv values for the original cube  
	double cubtv[10];
	
	//initializing the vector cubtv and VECGRIDRECT
	for (i=1 ; i < 10; i++ ){
		
		cubtv[i] = 0.0;
	
	}

	for (j=0; j < 10 ; j++) {
	
		VECMGRIDRECT_A[j] = 0.0;
		VECMGRIDRECT_B[j] = 0.0;
			
	}
	
	//Initializing variables of the integration points number
	SUM=0.0;
	expt=0.0;
	intresult = 0.0;
	initotal = 0.0;
	freeze = 0;
	
	//initializing the vector SUMdelt
	for (i=1 ; i <= 4; i++ ){
		
		SUMdelt[i] = 0.0;
	
	}
	
	//Attribution of values for cubtv
	for (i=1 ; i < 4; i++ ){
		
		cubtv[i] = GRIDRECT3D[3 + i] - GRIDRECT3D[i];
	
	}
	
	volume = cubtv[1]*cubtv[2]*cubtv[3];
	volume = (volume/NT_MINICUB);

	//initializing the matrix SUMdelt
	for (i=0 ; i < 4; i++ ){
		for (j=0 ; j <= NFAT[1]; j++ ){

			FAT_INTGR[i][j]=0.0;
				
		}
	}
	
	//Allocating the vector of the volume
	VOLUME =  F_VEC_ALLOC(NT_MINICUB+1);

	//Initializing the volume vector
	for (i=0; i <= NT_MINICUB; i++){
		VOLUME[i] = 0.0;	
	}

	for (i = 1; i <= NT_MINICUB; i++){
	
		VOLUME[i] = (MINIRECTSIZE[i][4] - MINIRECTSIZE[i][1])*(MINIRECTSIZE[i][5] - MINIRECTSIZE[i][2])*(MINIRECTSIZE[i][6] - MINIRECTSIZE[i][3]);
		
		if (VOLUME[i] <= volume*(0.1) ){
			
			freeze = 1;
			
		}
		
	}
	expt=-2.0;

	//FOR AXIL X
	//---------------------------------------------------------------------------------------------------------
	for (a=1; a <= NFAT[1]; a++){
		for (b=1; b <= NFAT[2]; b++){
			for (c=1; c <= NFAT[3]; c++){
			
				FAT_INTGR[1][a] += VECINT3D[a][b][c];
				
				//printf("\n VECINT: %15.6lf  FAT_INTGR: %15.6lf", VECINT3D[a][b][c], FAT_INTGR[1][a]);
				
			} 
		}

		//if(FAT_INTGR[1][a] <= (fatnorm) ){		

		//	FAT_INTGR[1][a] = (fatnorm);

		//}
	}
	//---------------------------------------------------------------------------------------------------------
	
	for (i=1; i<= NT_MINICUB; i++){
		
		SUM += VECINT[i];
		
	}

	//Calculating the new deltx
	for (i=1; i <= NFAT[1]; i++){
				
		FAT_INTGR[1][i] = FAT_INTGR[1][i]/SUM;	

	}

	//Calculating the new deltx
	for (i=1; i <= NFAT[1]; i++){
				
		MINICUBSIZE[1][i] = MINICUBSIZE[1][i]*(exp (expt*fabs(dbreak)*fabs(FAT_INTGR[1][i]) ));
		SUMdelt[1] += MINICUBSIZE[1][i];

	}
	
	//Normalizing the new deltx
	for (j=1; j <= NFAT[1]; j++){
				
		MINICUBSIZE[1][j] = (cubtv[1] * MINICUBSIZE[1][j])/SUMdelt[1];

	}
	
	//Generating the coordinates of matrix for each mini-cub
	//Attribution the coordinates of the slicing in the axil x
	for (i=1; i <= NFAT[1]; i++){
			
		FAT_EAXIL[1][i] = FAT_EAXIL[1][i-1] + MINICUBSIZE[1][i];

	}
	//Attribution the coordinates of the slicing in the axil y
	for (i=1; i <= NFAT[2]; i++){
			
		FAT_EAXIL[2][i] = FAT_EAXIL[2][i-1] + MINICUBSIZE[2][i];

	}
	//Attribution the coordinates of the slicing in the axil z
	for (i=1; i <= NFAT[3]; i++){
			
		FAT_EAXIL[3][i] = FAT_EAXIL[3][i-1] + MINICUBSIZE[3][i];

	}
	
	if(freeze == 0){
		//Initializing the counter
		count=1;
		
		//Generating the coordinates of each minicub
		for (i=0; i < NFAT[1]; i++){
			for (j=0; j < NFAT[2]; j++){
				for (n=0; n < NFAT[3]; n++){
					
					MINIRECTSIZE[count][1]= FAT_EAXIL[1][i];
					MINIRECTSIZE[count][2]= FAT_EAXIL[2][j];
					MINIRECTSIZE[count][3]= FAT_EAXIL[3][n];
					MINIRECTSIZE[count][4]= FAT_EAXIL[1][i+1];
					MINIRECTSIZE[count][5]= FAT_EAXIL[2][j+1];
					MINIRECTSIZE[count][6]= FAT_EAXIL[3][n+1];
					
					count++;		
						
				}
			}	
		}
	}

	//FOR AXIL Y
	//---------------------------------------------------------------------------------------------------------
	for (b=1; b <= NFAT[2]; b++){  
		for (a=1; a <= NFAT[1]; a++){   
			for (c=1; c <= NFAT[3]; c++){   
			
				FAT_INTGR[2][b] += VECINT3D[a][b][c] ;
	
			
			} 
		}

	}
	//---------------------------------------------------------------------------------------------------------
	
	for (i=1; i<= NT_MINICUB; i++){
		
		SUM += VECINT[i];
		
	}

	//Calculating the new deltx
	for (i=1; i <= NFAT[2]; i++){
				
		FAT_INTGR[2][i] = FAT_INTGR[2][i]/SUM;	

	}

	//Calculating the new deltx
	for (i=1; i <= NFAT[1]; i++){
				
		MINICUBSIZE[2][i] = MINICUBSIZE[2][i]*(exp (expt*fabs(dbreak)*fabs(FAT_INTGR[2][i]) ));
		SUMdelt[2] += MINICUBSIZE[2][i];

	}
	
	//Normalizing the new delty
	for (j=1; j <= NFAT[2]; j++){
				
		MINICUBSIZE[2][j] = (cubtv[2] * MINICUBSIZE[2][j])/SUMdelt[2];

	}
	
	//Generating the coordinates of matrix for each mini-cub
	//Attribution the coordinates of the slicing in the axil x
	for (i=1; i <= NFAT[1]; i++){
			
		FAT_EAXIL[1][i] = FAT_EAXIL[1][i-1] + MINICUBSIZE[1][i];

	}
	//Attribution the coordinates of the slicing in the axil y
	for (i=1; i <= NFAT[2]; i++){
			
		FAT_EAXIL[2][i] = FAT_EAXIL[2][i-1] + MINICUBSIZE[2][i];

	}
	//Attribution the coordinates of the slicing in the axil z
	for (i=1; i <= NFAT[3]; i++){
			
		FAT_EAXIL[3][i] = FAT_EAXIL[3][i-1] + MINICUBSIZE[3][i];

	}
	
	//Initializing the counter
	count=1;
	intresult = 0.0;
	
	if(freeze == 0){
		//Generating the coordinates of each minicub
		for (i=0; i < NFAT[1]; i++){
			for (j=0; j < NFAT[2]; j++){
				for (n=0; n < NFAT[3]; n++){
					
					MINIRECTSIZE[count][1]= FAT_EAXIL[1][i];
					MINIRECTSIZE[count][2]= FAT_EAXIL[2][j];
					MINIRECTSIZE[count][3]= FAT_EAXIL[3][n];
					MINIRECTSIZE[count][4]= FAT_EAXIL[1][i+1];
					MINIRECTSIZE[count][5]= FAT_EAXIL[2][j+1];
					MINIRECTSIZE[count][6]= FAT_EAXIL[3][n+1];
					
					count++;		
						
				}
			}	
		}	
	}
	
	//FOR AXIL Z
	//---------------------------------------------------------------------------------------------------------
	for (c=1; c <= NFAT[3]; c++){
		for (a=1; a <= NFAT[1]; a++){
			for (b=1; b <= NFAT[2]; b++){
			
				FAT_INTGR[3][c] += VECINT3D[a][b][c] ;
						
			} 
		}

	}
	//---------------------------------------------------------------------------------------------------------

	for (i=1; i<= NT_MINICUB; i++){
		
		SUM += VECINT[i];
		
	}

	//Calculating the new deltx
	for (i=1; i <= NFAT[3]; i++){
				
		FAT_INTGR[3][i] = FAT_INTGR[3][i]/SUM;	

	}

	//Calculating the new deltz
	for (i=1; i <= NFAT[3]; i++){
				
		MINICUBSIZE[3][i] = MINICUBSIZE[3][i]*(exp (expt*fabs(dbreak)*fabs(FAT_INTGR[3][i]) ));
		SUMdelt[3] += MINICUBSIZE[3][i];

	}
	
	//Normalizing the new deltz
	for (j=1; j <= NFAT[3]; j++){
				
		MINICUBSIZE[3][j] = (cubtv[3] * MINICUBSIZE[3][j])/SUMdelt[3];

	}
	
	//Generating the coordinates of matrix for each mini-cub
	//Attribution the coordinates of the slicing in the axil x
	for (i=1; i <= NFAT[1]; i++){
			
		FAT_EAXIL[1][i] = FAT_EAXIL[1][i-1] + MINICUBSIZE[1][i];

	}
	//Attribution the coordinates of the slicing in the axil y
	for (i=1; i <= NFAT[2]; i++){
			
		FAT_EAXIL[2][i] = FAT_EAXIL[2][i-1] + MINICUBSIZE[2][i];

	}
	//Attribution the coordinates of the slicing in the axil z
	for (i=1; i <= NFAT[3]; i++){
			
		FAT_EAXIL[3][i] = FAT_EAXIL[3][i-1] + MINICUBSIZE[3][i];

	}
	
	//Initializing the counter
	count=1;
	
	if(freeze == 0){
		//Generating the coordinates of each minicub
		for (i=0; i < NFAT[1]; i++){
			for (j=0; j < NFAT[2]; j++){
				for (n=0; n < NFAT[3]; n++){
					
					MINIRECTSIZE[count][1]= FAT_EAXIL[1][i];
					MINIRECTSIZE[count][2]= FAT_EAXIL[2][j];
					MINIRECTSIZE[count][3]= FAT_EAXIL[3][n];
					MINIRECTSIZE[count][4]= FAT_EAXIL[1][i+1];
					MINIRECTSIZE[count][5]= FAT_EAXIL[2][j+1];
					MINIRECTSIZE[count][6]= FAT_EAXIL[3][n+1];
					
					count++;		
						
				}
			}	
		}	
	}
	
	//printf("\n\n\n        -------------------------------------------------------------------------");
	//printf("\n        ------------       GENERATING MINIRECTSIZE PARAMETRS       --------------");
	//printf("\n        -------------------------------------------------------------------------\n\n\n");
	//printf("        Coordinates of the points that define of the each AFTER ADAPTATION mini-cubs:\n");
	//printf("        -------------------------------------------------------\n\n");
	//printf("\n\n1");
	//for (i = 1; i <= NT_MINICUB; i++) {
	//	printf("        X    %15.6lf %15.6lf %15.6lf \n",
	//			MINIRECTSIZE[i][1], MINIRECTSIZE[i][2], MINIRECTSIZE[i][3] ); 
	//	printf("        X    %15.6lf %15.6lf %15.6lf \n",
	//			MINIRECTSIZE[i][4], MINIRECTSIZE[i][5], MINIRECTSIZE[i][6] );
				
	
	
	F_VEC_DEALLOC(VOLUME);
}


#endif




