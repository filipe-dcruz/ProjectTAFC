#ifndef INPUT_HEADER
#define INPUT_HEADER

#define NSPE 4 /*Number of species*/

struct Species{
  const char* name ;
  double rqm ;
} ;

struct Species Specie[NSPE] ;

//For simplicity the initial values for the box and time are zero

// particles parameters
static const int NPAR = 10 ;

// Spatial parameters
#define NX 100
static const double BOX = 1.0 ;

// Time parameters
static const double TMAX = 1.0 ;
static const double CVAL = 1.0 ; // Courant number

// Diagnostics
#define NDUMP 50

// Declaration of initial values
static inline double BxInicial( double x ){ return 0.; }
static inline double ByInicial( double x ){ return 0.; }
static inline double BzInicial( double x ){
  return 0.;
}

static inline double ExInicial( double x ){
  return 0. ;
}
static inline double EyInicial( double x ){ return 0.; }
static inline double EzInicial( double x ){ return 0.; }

#endif
