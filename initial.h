/**
    Declaration and definition of multiple calculated parameters and functions
    responsable for the initiation of the variables and objects.
    @file initial.h
    @author Filipe Cruz
**/
#ifndef INITIAL_HEADER
#define INITIAL_HEADER

#include <iostream>
#include <random>
#include <cstdlib>
#include <ctime>
#include <cmath>

#include <sys/stat.h>
#include <sys/types.h>

#include "input.h"

#define NPARLIMIT 1000      /* Limit for the number of particle per cell*/
#define NFIELDS 6           /*Aux : number of electromagnetic fields*/
#define MKDIR 6             /*Aux : length of mkdir command*/
#define RM1 6               /*Aux : length of rm command*/
#define RM2 5               /*Aux : length of rm command*/

#define MAX_FILE_NAME 30    // Maximum size for the output files name

// Name of directory
extern char* name_file ;    // Name of the directory for the output files

// Simulation parameters
static const double BOX = XF-X0 ;        // Size of the simulation box
static const double dx = BOX/NX ;        // Width of cells
static const double dt = dx*CVAL ;       // Time step

// Grids for the electromagnetic fields
extern double* Ex[NDIM], *Bx[NDIM] ;
// Positions of the grid cells and their borders
extern double* xx ;

// Name of the output file with the simulation parameter values
extern char output_file[MAX_FILE_NAME] ;

// Name of the output files for the fields components
extern char field_files[NFIELDS][MAX_FILE_NAME];
// Pointers to the correspondent eletromagnetic fields
extern double* field_var[NFIELDS] ;

// Name of the output files for the species density
extern char density_files[NSPE][MAX_FILE_NAME] ;
// Attachment for the the name of these files
static const char density_name[MAX_FILE_NAME] = "-density.txt" ;
// Name of the output files for the species phase space
extern char val_files[NSPE][MAX_FILE_NAME] ;
// Attachment for the the name of these files
static const char val_name[MAX_FILE_NAME] = "-val.txt" ;

// Auxiliary parametes - calculated for efficiency
static const int NX1 = NX-1 , NX2 = NX-2, NXp1 = NX+1 ;
static const double dif1 = 2*M_PI/NX , dif2 = dif1/2. ;
static const double aux_ = 2*M_PI/BOX , aux1_ = 1/dx ;
//, _dx2 = 1/(dx*2.) ;

// Auxiliary Functions
// Check if the parameter of the simulation are valid quantities.
int CheckParameters( void ) ;

// Create names for output files
void InitialDefinitions( const char * ) ;

// Calculate the initial values for the grids and species
void DefineInitialValues( void ) ;

// Calculate the density from the particles position
void CalculateTheDensity( void ) ;

// Release temporary allocated memory of the grids
void FinalDeclarations( void ) ;

#endif
