#include "species.h"

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

}
