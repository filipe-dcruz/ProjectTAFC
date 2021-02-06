#ifndef INITIAL_HEADER
#define INITIAL_HEADER

#include <stdio.h>
#include <stdlib.h> /*exit*/

# include "input.h"


// Secundaty parameters
double dx , dt ; // Spatial and time steps

//Declare array for particles
double * xpar , *vpar ;

//Declare grid
double *Ex, *Ey, *Ez ;
double *Bx, *By, *Bz ;

// Auxiliary functions
void InitialDeclarations( void ) ;
void DefineInitialValues( void ) ;

void FinalDeclarations( void ) ;

#endif
