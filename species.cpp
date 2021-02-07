#include "species.h"

// Auxiliary functions
void Species::CalculateLastParameters( double dx )
{
  NParTot = int((xf-x0)/dx)*npar ; //Total number of particles
  qc = den/npar ;
}

/***Definitions of Species methods***/
//Constuct
Species::Species(const char *_name ,double _rqm, uint _npar ,
  double _v0, double _vth,double _x0, double _xf, double _den ):
  rqm(_rqm) , npar(_npar) , v0(_v0), vth(_vth) , x0(_x0),
  xf(_xf), den(_den)
{
  for ( int i = 0 ; i < NAME_LIMIT && _name[i] ; i++ ){
    name[i] = _name[i] ;
  }
}

//Destructor
Species::~Species()
{
  delete xval ;
  delete pval ;

  //free(xpos) ;
  //free(ppos) ;

}

void Species::CreateList( uint N , double dx )
{

  // Calculate remaining parameters
  CalculateLastParameters(dx) ;

  //size_t aux   = N * sizeof(itr*) ;

  // Create list of particles
  //xpos = ( itr* ) malloc(aux) ;
  xval = new std::list<double> ( NParTot ) ;

  // Create list of velocities
  //ppos = ( itr* ) malloc(aux) ;
  pval = new std::list<double> ( NParTot ) ;

}
