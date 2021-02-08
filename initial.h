#ifndef INITIAL_HEADER
#define INITIAL_HEADER

#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>

#include <sys/stat.h>
#include <sys/types.h>

#include "input.h"

#define NPARLIMIT 1000
#define NFIELDS 6
#define MKDIR 6
#define RM1 6
#define RM2 5

#define MAX_FILE_NAME 30

// Auxiliary
const int NX1 = NX-1 ;

// Name of directory
extern char* name_file ;

// Secundaty parameters
static const double BOX = XF-X0 ;
static const double dx = BOX/NX ;
static const double dt = dx*CVAL ; // Spatial and time steps

//Declare grid for electromagnetic fields
extern double* Ex[NDIM], *Bx[NDIM] ;
extern double* xx ;

extern char output_file[MAX_FILE_NAME] ;

extern char field_files[NFIELDS][MAX_FILE_NAME];
extern double* field_var[NFIELDS] ;

extern char density_files[NSPE][MAX_FILE_NAME] ;
static const char density_name[MAX_FILE_NAME] = "-density.txt" ;
extern char val_files[NSPE][MAX_FILE_NAME] ;
static const char val_name[MAX_FILE_NAME] = "-val.txt" ;

// Auxiliary Functions
int CheckParameters( void ) ;
void InitialDefinitions( const char * ) ;
void DefineInitialValues( void ) ;

void FinalDeclarations( void ) ;

#endif
