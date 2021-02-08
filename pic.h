#ifndef PIC_HEADER
#define PIC_HEADER

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "input.h"
#include "initial.h"

// Calculate cross product between two vectors
inline void CrossProduct(double* v1, double* v2, double* vr){
  vr[0] = v1[1]*v2[2]-v1[2]*v2[1] ;
  vr[1] = v1[2]*v2[0]-v1[0]*v2[2] ;
  vr[2] = v1[0]*v2[1]-v1[1]*v2[0] ;
}

bool PrintDiagnostics(double,double***) ;

void CalculateNewPosVel( double*** , double*** ) ;

void UpdateData( double*** , double*** ) ;

void ComputePIC( const char *) ;

void ProduceDiagnostics(void) ;

void PIC1D( const char *) ;

#endif
