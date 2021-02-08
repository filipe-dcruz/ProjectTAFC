#ifndef INPUT_HEADER
#define INPUT_HEADER

#include "species.h"

#define NSPE 4  /*Number of species*/
#define NCOR 4  /*Number of cores  */

static const double CVAL = 0.995 ; // Courant number

// Spatial parameters
#define NX 600                      // Number of cells
static const double X0 = 0.0 ;
static const double XF = 6.0 ;

// Time parameters
static const double TMIN = 0.0 ;
static const double TMAX = 100.0 ;

// Diagnostics
const uint NDUMP = 100 ; //Interactions
const double DUMP_PER = 0.05 ;

// Declaration of initial values
namespace InitialFields
{
  double BxInicial( double x ) ;
  double ByInicial( double x ) ;
  double BzInicial( double x ) ;

  double ExInicial( double x ) ;
  double EyInicial( double x ) ;
  double EzInicial( double x ) ;
}

extern Species specie[NSPE] ;

#endif
