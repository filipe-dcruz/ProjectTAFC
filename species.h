#ifndef SPECIES_HEADER
#define SPECIES_HEADER

#include <list>
#include <cstdlib>
#include <iostream>

#define NAME_LIMIT 20
#define NDIM 3

typedef unsigned int uint ;
typedef std::list<double>::iterator itr ;
typedef std::list<double> type ;

/*
  Define Species that contains the information of the initial and contant
  values for each specie.
*/
class Species{
  uint NParTot ;
  double qc ;
  double ql ;

  void CalculateLastParameters ( double , double , double ) ;

public:
  // Variables
  // Global
  char name [NAME_LIMIT]; // Name of species
  const double rqm ;       // q/m ratio
  const uint npar ; //Number of particles per cell of this specie

  // Velocity
  const double v0[NDIM] ;
  const double vth[NDIM] ;

  double x0, xf ;
  const double den ;

  // Arrays
  std::list<double>* pval[NDIM] ;
  std::list<double>* xval[NDIM] ;

  double * density ;

  //Functions
  Species( const char*, const double, const uint, const double*, const double*,
    const double, const double, const double ) ;
  ~Species () ;

  void CreateList( uint , double , double , double ) ;

  uint NumOfPar(){ return NParTot ; };
  double qlvalue(){ return ql ; };
  double qcvalue(){ return qc ; };
};

#endif
