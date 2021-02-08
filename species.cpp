#include "species.h"

// Auxiliary functions
void Species::CalculateLastParameters( double dx , double dt )
{
  NParTot = int((xf-x0)/dx)*npar ; //Total number of particles
  qc = den/npar ;
  ql = rqm/2.*dt ;
}

/***Definitions of Species methods***/
//Constuct
Species::Species(const char *_name , const double _rqm, const uint _npar ,
  const double* _v0, const double* _vth,
  const double _x0, const double _xf, const double _den ):
  rqm(_rqm), npar(_npar),
  v0{_v0[0],_v0[1],_v0[2]}, vth{_vth[0],_vth[1],_vth[2]},
  x0(_x0), xf(_xf), den(_den)
{
  for ( int i = 0 ; i < NAME_LIMIT && _name[i] ; i++ )
    name[i] = _name[i] ;
}

//Destructor
Species::~Species()
{
  for ( int i = 0 ; i < NDIM ; i++ ){
    delete xval[i] ;
    delete pval[i] ;
  }

  free(density) ;

}

void Species::CreateList( uint N , double dx , double dt )
{
  // Calculate remaining parameters
  CalculateLastParameters(dx,dt) ;

  // Create list of particles and velocities
  for( int i = 0 ; i < NDIM ; i++ ){
    xval[i] = new std::list<double> ( NParTot ) ;
    pval[i] = new std::list<double> ( NParTot ) ;
  }

  // Calculate the density
  density = new double[N] ;
}
