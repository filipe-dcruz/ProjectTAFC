#include "initial.h"

/*
  Initiate parameters for the simulation ;
*/
void InitialDeclarations( ){

	// Check invalid values for parameters
  if ( !( NPAR > 0 ) ){
		printf("--ERROR: INVALID NUMBER OF PARAMETERS\n") ;
		exit(1);
	}
  if ( !( NX > 0 ) ){
		printf("--ERROR: INVALID GRID NUMBER\n") ;
		exit(1);
	}
  if ( !( BOX > 0. ) ){
		printf("--ERROR: INVALID \"BOX\" VALUE\n") ;
		exit(1);
	}
  if ( !( TMAX > 0. ) ){
		printf("--ERROR: INVALID \"TMAX\" VALUE\n") ;
		exit(1);
	}
  if ( !( CVAL > 0. ) ){
		printf("--ERROR: INVALID COURANT PARAMETER\n") ;
		exit(1);
	}
	if ( CVAL > 1. )
		printf("--WARNING: Courant number is lower than 1.0\n") ;

  //Calculate steps and auxiliary functions#include "initial.h" //Input configuration
  dx = BOX / ( (double) NX ) ;
  dt = CVAL * dx ;

  // Define grids
  Ex = (double *) malloc( NX * sizeof(double));
  Bz = (double *) malloc( NX * sizeof(double));

  // Define particles
  xpar = (double *) malloc( NPAR * sizeof(double));
  vpar = (double *) malloc( NPAR * sizeof(double));
}

/*
  Define parameters with the setup of the simulation ;
*/
void DefineInitialValues(){
  double x = 0. ;

  // initiate grids
  for( int i = 0 ; i < NX ; i++ , x += dx ){
    Bz[i] = BzInicial( x ) ;
    Ex[i] = ExInicial( x ) ;
  }
}

void FinalDeclarations( ){
  printf("NX = %d\n", (int) NX) ;
}
