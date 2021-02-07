#ifndef SPECIES_HEADER
#define SPECIES_HEADER

#include <list>
#include <cstdlib>
#include <iostream>

#define NAME_LIMIT 20

typedef unsigned int uint ;
typedef std::list<double>::iterator itr ;

/*
  Define Species that contains the information of the initial and contant
  values for each specie.
*/
class Species{
  uint NParTot ;
  double qc ;

  void CalculateLastParameters ( double ) ;

public:
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

  // Arrays
  std::list<double>* pval ;
  std::list<double>* xval ;

  itr** xpos ;
  itr** ppos ;

  //Functions
  Species ( const char*, double, uint, double, double,
    double, double, double ) ;
  ~Species () ;

  void CreateList( uint , double ) ;

  uint NumOfPar(){ return NParTot ; };
} ;

#endif
