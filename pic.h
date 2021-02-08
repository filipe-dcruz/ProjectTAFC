#ifndef PIC_HEADER
#define PIC_HEADER

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "input.h"
#include "initial.h"

bool PrintDiagnostics(double) ;

void CalculateNewPosVel() ;

void ComputePIC( const char *) ;

void ProduceDiagnostics(void) ;

void ComputePosVel(void) ;

void PIC1D( const char *) ;

#endif
