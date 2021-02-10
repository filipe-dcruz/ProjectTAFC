/**
    Declares the class Species that contains the information for each specie
    of particles.
    This declaration was not necessary, but it was used to better organize the
    code and to use the private attributes.

    @file species.h
    @author Filipe Cruz
*/
#ifndef SPECIES_HEADER
#define SPECIES_HEADER

#include <cstdlib>
#include <iostream>

#define NAME_LIMIT 20       /*Maximum size for the name of the particles*/
#define NDIM 3              /*Number of dimensions-Auxiliar*/

//typedef unsigned int uint ;

/*
  Define Species that contains the information of the initial configurations
  and contant values for each specie od particles.
*/
class Species{
  int NParTot ;             // Total number of particles
  double qc ;               // Cloud charge
  double ql ;               // q' value for the Boris pusher

  // Calculates the values of auxiliar parameters
  void CalculateLastParameters ( double , double , double ) ;

public:
  // Variables
  // Global
  char name [NAME_LIMIT];   // Name of species
  const double rqm ;        // q/m ratio
  const int npar ;          // Number of particles per cell

  // Velocity
  const double v0[NDIM] ;   // Fluid velocities
  const double vth[NDIM] ;  // Thermal vecloities

  double x0, xf ;           // Vondaries of the uniform distibution
  const double den ;        // Initial uniform density

  // Arrays
  double* pval[NDIM] ;      // Mommentum for each particles - dynamic
  double* xval[NDIM] ;      // Position for each particles - dynamic

  double * density ;        // Density distribution - dynamic

  // Costructor for the species
  Species( const char*, const double, const int, const double*, const double*,
    const double, const double, const double ) ;
  ~Species () ;

  // Allocates the memory for the particles
  void CreateList( int , double , double , double ) ;

  // Auxiliary functions. They were created to avoid an acidental change of
  // these calculated parameters
  int NumOfPar(){ return NParTot ; };
  double qlvalue(){ return ql ; };
  double qcvalue(){ return qc ; };
};

#endif
