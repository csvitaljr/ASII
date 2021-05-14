//
//  asii_test.h
//
//  Created by Renaldo Tenorio de Moura Junior and Carlos Vital dos Santos Junior on 13/05/21.
//  Copyright © 2021 Renaldo Tenorio de Moura Junior. All rights reserved.
//
//  This function adapt the  coordinates that define the size of each mini-cub for a good parallel balanncig

#ifndef asii_test_h
#define asii_test_h

double ASII_TEST(double CUR_POINT[4]){

    double R2;
    double X=0.0;
    double Y=0.0;
    double Z=0.0;
    double N=(1/39.6876);
    double a=-0.27001;
    double gauss;
    
    X=CUR_POINT[1];
    Y=CUR_POINT[2];
    Z=CUR_POINT[3];

    R2=(X*X)+(Y*Y)+(Z*Z);
    gauss=N*exp(a*R2);
    //gauss=1.0;

    return(gauss);

}


#endif /* asii_test.h */
