#ifndef INITIAL_HEADER
#define INITIAL_HEADER

#include <iostream>
#include <random>

#include "input.h"

#define NPARLIMIT 1000

// Secundaty parameters
static const double BOX = XF-X0 ;
static const double dx = BOX/NX ;
static const double dt = dx*CVAL ; // Spatial and time steps

//Declare grid for electromagnetic fields
extern double *Ex, *Ey, *Ez ;
extern double *Bx, *By, *Bz ;
extern double *xx ;

// Auxiliary Functions
int CheckParameters() ;
void InitialDefinitions( void ) ;
void DefineInitialValues( void ) ;

void FinalDeclarations( void ) ;

#endif
