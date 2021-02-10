/**
    Defines the class Species that contains the information for each specie
    of particles. Declarations made in file specs.h

    @file species.cpp
    @author Filipe Cruz
*/

#include "species.h"

/*
  Calculates constnat parameters that will be used. It increases efficiency
  @param
    dx  : distance between cells
    dt  : time step
    X0_ : Left boundary spatial value of the simulation box
*/
void Species::CalculateLastParameters( double dx , double dt , double X0_ )
{
  // Total number of particles
  NParTot = int((xf-x0)/dx)*npar ; // Number of cells * npar

  // Calculate cloud charge of the particle.
  // Negative for negative charged particles
  if ( rqm > 0 )
    qc = den/npar ;
  else
    qc = -den/npar ;

  // Auxiliar paramenter for the Boris pusher
  ql = rqm/2.*dt ;

  // Subtraxts the left boundary of the particle to increase efficiency
  xf -= X0_ ;
  x0 -= X0_ ;
}

/*
  Contruct for the Species class
  @param
    _name   : name of the specie
    _rqm    : ratio charge/mass
    _npar   : number of particles per cell
    _v0     : fluid velocity in three directions
    _vth    : thermal velocity in three directions
    _x0     : initial left spatial value for the specie distibution
    _xf     : initial right spatial value for the specie distibution
    _den    : density of the uniform distribution
*/
Species::Species(const char *_name , const double _rqm, const uint _npar ,
  const double* _v0, const double* _vth,
  const double _x0, const double _xf, const double _den ):
  rqm(_rqm), npar(_npar),
  v0{_v0[0],_v0[1],_v0[2]}, vth{_vth[0],_vth[1],_vth[2]},
  x0(_x0), xf(_xf), den(_den)
{
  // Sets the species name
  for ( int i = 0 ; i < NAME_LIMIT && _name[i] ; i++ )
    name[i] = _name[i] ;
}

/*
  Destructor of the Species class. Frees allocated memory
*/
Species::~Species()
{
  // Free phase space memory
  for ( int i = 0 ; i < NDIM ; i++ ){
    delete xval[i] ; delete pval[i] ;
  }

  // Free density memory
  free(density) ;
}


/*
  Defines the parameters dependent on the simulation paramenters and allocates
  the necessary memory
  @param
    N    : number of cells
    dx  : distance between cells
    dt  : time step
    X0_ : Left boundary spatial value of the simulation box
*/
void Species::CreateList( int N , double dx , double dt , double X0_ )
{

  // Calculate remaining and auxiliar parameters for the specie
  CalculateLastParameters(dx,dt,X0_) ;

  // Allocates memory for the position and velocities of the particles
  for( int i = 0 ; i < NDIM ; i++ ){
    xval[i] = new double [NParTot] ;
    pval[i] = new double [NParTot] ;
  }

  // Allocates memory for the density of the particles
  density = new double[N] ;
}
