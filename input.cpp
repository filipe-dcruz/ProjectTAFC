/**
    Definition of the paramenters for the configuration of the PIC simulation.
    @file input.cpp
    @author Filipe Cruz
**/

#include "input.h"

// x component function of the initial magnetic field
double InitialFields::BxInicial( double x ) { return 0. ; }

// y component function of the initial magnetic field
double InitialFields::ByInicial( double x ) { return 0. ; }

// z component function of the initial magnetic field
double InitialFields::BzInicial( double x )
{
  if ( x < 3.0 )
    return 0. ;
  else
    return 0.0 ;
}

// x component function of the initial electric field
double InitialFields::ExInicial( double x ) { return 0. ; }

// y component function of the initial electric field
double InitialFields::EyInicial( double x ) { return 0. ; }

// z component function of the initial electric field
double InitialFields::EzInicial( double x ){ return 0. ; }

//Auxiliary vectors for the initiation of the velocities
static const double v0Aux1[NDIM] = {0.1,0.,0.} ;
static const double v0Aux2[NDIM] = {-0.1,0.,0.} ;
static const double vthAux1[NDIM] = {0.001,0.001,0.001} ;
static const double vthAux2[NDIM] = {0.01,0.01,0.01} ;

// Declares and initiates the species for the simulation 
Species specie[NSPE] = {
   Species("ions"        ,100.,10,
    v0Aux1,vthAux1,5.0,6.0,1.),
   Species("electrons"   ,-1. ,10,
    v0Aux1,vthAux2,5.0,6.0,1.),
   Species("bg-ions"     ,100.,10,
    v0Aux2,vthAux1,3.0,6.0,1.),
   Species("bg-electrons",-1. ,10,
    v0Aux2,vthAux2,3.0,6.0,1.)
};
