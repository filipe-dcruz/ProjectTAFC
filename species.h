#ifndef SPECIES_HEADER
#define SPECIES_HEADER

#define NAME_LIMIT 20

typedef unsigned int uint ;

/*
  Define Species that contains the information of the initial and contant
  values for each specie.
*/
struct Species{
  // Variables
  // Global
  char name [NAME_LIMIT]; // Name of species
  const double rqm ;       // q/m ratio
  const uint npar ; //Number of particles per cell of this specie

  // Velocity
  const double v0 ;
  const double vth ;

  const double x0, xf ;
  const double den ;

  // Arrayes
  double* pval ;
  double* xval ;

  //Functions
  Species ( const char*, double, uint, double, double,
    double, double, double ) ;
  ~Species () ;
} ;

#endif
