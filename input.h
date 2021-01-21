#ifndef INPUT_HEADER
#define INPUT_HEADER

//For simplicity the initial values for the box and time are zero

// particles parameters
static const int NPAR = 10 ;
static const int NSPE = 4 ;

// Spatial parameters
static const int NX = 100 ;
static const double BOX = 1.0 ;

// Time parameters
static const double TMAX = 1.0 ;
static const double CVAL = 1.0 ; // Courant number

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
