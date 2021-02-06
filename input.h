#ifndef INPUT_HEADER
#define INPUT_HEADER

#include "species.h"

#define NSPE 4  /*Number of species*/
#define NCOR 4  /*Number of cores  */

static const double CVAL = 1.0 ; // Courant number

// Spatial parameters
#define NX 100                      // Number of cells
static const double X0 = 0.0 ;
static const double XF = 1.0 ;

// Time parameters
static const double TMIN = 0.0 ;
static const double TMAX = 1.0 ;

// Diagnostics
const uint NDUMP = 50 ; //Interactions

// Declaration of initial values
namespace InitialFields
{
  static inline double BxInicial( double x ) ;
  static inline double ByInicial( double x ) ;
  static inline double BzInicial( double x ) ;

  static inline double ExInicial( double x ) ;
  static inline double EyInicial( double x ) ;
  static inline double EzInicial( double x ) ;
}

extern Species specie[NSPE] ;

#endif
