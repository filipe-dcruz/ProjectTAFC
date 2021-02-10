/**
    Files responsible by the particle in cell algorith
    @file pic.h
    @author Filipe Cruz
**/

#ifndef PIC_HEADER
#define PIC_HEADER

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "input.h"
#include "initial.h"

// Calculate cross product between two vectors
inline void CrossProduct(double* v1, double* v2, double* vr){
  vr[0] = v1[1]*v2[2]-v1[2]*v2[1] ;
  vr[1] = v1[2]*v2[0]-v1[0]*v2[2] ;
  vr[2] = v1[0]*v2[1]-v1[1]*v2[0] ;
}

// Print the results for a time step in the output files
bool PrintDiagnostics(double,double***) ;

// Calculates the phase space data for the next iteration
void CalculateNewPosVel( double*** , double*** ) ;

// Update the phase data with the next iteration data
void UpdateData( double*** , double*** ) ;

// Get vectors that will be used to calculate the field
void GetFourierVectors( double* ) ;

// Calculates the next values for the electric field
void FieldSolver( double * , double * ) ;

// Compute Particle-In-Cell the particle in cell algorith
void ComputePIC( const char *) ;

// Function responsable in createing the file with the configuration parameters
void ProduceDiagnostics(void) ;

// Function responsable in performing the 1D electrostatic PIC code
void PIC1D( const char *) ;

#endif
