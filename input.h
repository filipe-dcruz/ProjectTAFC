/**
    Declares the paramenters for the configuration of the PIC simulation.
    @file input.h
    @author Filipe Cruz
**/
#ifndef INPUT_HEADER
#define INPUT_HEADER

#include "species.h"

#define NSPE 2                          /*Number of species*/
#define NCOR 4                          /*Number of cores - NOT IMPLEMENTED */

static const double CVAL = 0.995 ;      // Default Courant number

// Spatial parameters
#define NX 600                          // Number of cells
static const double X0 = 0.0 ;          // Left boundary
static const double XF = 6.0 ;          // Right boundary

// Time parameters
static const double TMIN = 0.0 ;        // Initial time
static const double TMAX = 50.0 ;       // Final time

// Diagnostics
const uint NDUMP = 5 ;                  // Frequency of the dump
const double DUMP_PER = 0.1 ;           // Percentage of particles to dump

// Declaration of functions for initial eletromagnetic fields
namespace InitialFields
{
  // Magnetic field components
  double BxInicial( double x ) ;
  double ByInicial( double x ) ;
  double BzInicial( double x ) ;

  // Electric field components
  double ExInicial( double x ) ;
  double EyInicial( double x ) ;
  double EzInicial( double x ) ;
}

// Declaration of the Species object
extern Species specie[NSPE] ;

#endif
