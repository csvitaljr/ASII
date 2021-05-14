//
//  asii_slice.h
//
//  Created by Renaldo Tenório de Moura Júnior and Carlos Vital dos Santos Júnior on 23/08/2018
//  Copyright © 2018 Renaldo Tenorio de Moura Junior. All rights reserved.
//
//  This function generate and do the calculate of the points coordinates that define the size of each mini-cub

#ifndef asii_slice_h
#define asii_slice_h

// NFAT is the number of the mini-cub that will be generate
// GRIDRECT3D contains the cub-coordinates

 void ASII_SLICE(int NFAT[4], double GRIDRECT3D[20], double **MINIRECTSIZE, double **MINICUBSIZE , int NT_MINICUBS, double CUBTV[4], double **FAT_EAXIL) {
		
	//Defining variables here used
	int i;
	int j;
	int n;
	int count;
	
	double deltxyz[4];
	double SUMdelt;
	
	//printf("0 Ok Until Here !! :) \n");
	
	// First put zero in all spaces. Assuming **MINIRECTSIZE has size [NFAT][7]
	for (i= 1 ; i < NT_MINICUBS +1 ; i++){
		for (n=0 ; n < 7 ; n++){
			
			MINIRECTSIZE[i][n] = 0.0;
		}
	}
	
	//Calculing the first delts
	deltxyz[1]= ( CUBTV[1] )/( NFAT[1] );
	deltxyz[2]= ( CUBTV[2] )/( NFAT[2] );
	deltxyz[3]= ( CUBTV[3] )/( NFAT[3] );	
	
	
	//printf("1 Ok Until Here !! :) \n");
	
	
	//Passing the values of delta x, y and z for the mini cubs size matrix 
	for (i=1; i <= NFAT[1]; i++){
		
		MINICUBSIZE[1][i] = deltxyz[1];

	}
	
	for (i=1; i <= NFAT[2]; i++){
		
		MINICUBSIZE[2][i] = deltxyz[2];

	}
	
	for (i=1; i <= NFAT[3]; i++){
		
		MINICUBSIZE[3][i] = deltxyz[3];

	}
	
	//printf("2 Ok Until Here !! :) \n");
	
	//Generating the coordinates of matrix for each mini-cub

	//Atribuition the point original of the cube 
	for (i=1; i <= 3; i++){
		
		FAT_EAXIL[i][0] = GRIDRECT3D[i];
		
	} 
	
	//Atribuition the coordinates of the slicing in the axil x
	for (i=1; i <= NFAT[1]; i++){
			
		FAT_EAXIL[1][i] = FAT_EAXIL[1][i-1] + MINICUBSIZE[1][i];

	}
	//Atribuition the coordinates of the slicing in the axil y
	for (i=1; i <= NFAT[2]; i++){
			
		FAT_EAXIL[2][i] = FAT_EAXIL[2][i-1] + MINICUBSIZE[2][i];

	}
	//Atribuition the coordinates of the slicing in the axil z
	for (i=1; i <= NFAT[3]; i++){
			
		FAT_EAXIL[3][i] = FAT_EAXIL[3][i-1] + MINICUBSIZE[3][i];

	}
	
	//printf("3 Ok Until Here !! :) \n");
	
	//Initializing the counter
	count=1;
	
	//Genarating the coordinates of each minicub
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
	
	//printf("4 Ok Until Here !! :) \n");
	
	NT_MINICUBS = NFAT[1]*NFAT[2]*NFAT[3];
	/*
	printf("\n\n\n        -------------------------------------------------------------------------");
	printf("\n        ------------       GENERATING MINIRECTSIZE PARAMETRS       --------------");
	printf("\n        -------------------------------------------------------------------------\n\n\n");
		
	
	printf("        The GRIDRECT Coordinates:\n");
	printf("        --------------------------\n\n");
	printf( "        BEGIN    %15.6lf %15.6lf %15.6lf \n",
				GRIDRECT3D[1], GRIDRECT3D[2], GRIDRECT3D[3] ); 
	printf( "        FINAL    %15.6lf %15.6lf %15.6lf \n\n",
				GRIDRECT3D[4], GRIDRECT3D[5], GRIDRECT3D[6] );
				
	printf("        The GRIDRECT tv's coordinates:\n");
	printf("        ------------------------------\n\n");	
	printf( "        Tv     %15.6lf           0.000000        0.000000 \n", CUBTV[1]);
	printf( "        Tv            0.000000    %15.6lf        0.000000 \n", CUBTV[2]);
	printf( "        Tv            0.000000           0.000000 %15.6lf \n", CUBTV[3]);
	printf(	"\n\n");
				
	printf("        Coordinates of the points that define of the mini-cubs:\n");
	printf("        -------------------------------------------------------\n\n");

		for (i = 1; i <= NT_MINICUBS; i++) {
	printf("        X    %15.6lf %15.6lf %15.6lf \n",
			MINIRECTSIZE[i][1], MINIRECTSIZE[i][2], MINIRECTSIZE[i][3] ); 

	printf("        X    %15.6lf %15.6lf %15.6lf \n",
			MINIRECTSIZE[i][4], MINIRECTSIZE[i][5], MINIRECTSIZE[i][6] );
				
	//fprintf(output, "\n");
	
	}
	
	printf("\n         Values of the Delts and Slicing:\n");
	printf("        -------------------------------------------------------");		
	
	for (i = 0; i < NFAT[1]+1; i++) {
		printf("\n");
		for(j=1; j <= 3; j++){
			
			printf(             "%15.6f %15.6f \n", MINICUBSIZE[j][i], FAT_EAXIL[j][i]);
		}
	}
	
	printf("\n");
				
	*/
	
	//printf("5 Ok Until Here !! :) \n");

}
	
#endif





