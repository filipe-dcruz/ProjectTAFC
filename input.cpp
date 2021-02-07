#include "input.h"

//Components of the initial magnetic field
double InitialFields::BxInicial( double x )
{
  return 0. ;
}

double InitialFields::ByInicial( double x )
{
  return 0 ;
}

double InitialFields::BzInicial( double x )
{
  if (x < 0.5 )
    return 0. ;
  else
    return 1.0 ;
}

//Components of the initial electric field
double InitialFields::ExInicial( double x )
{
  return 0. ;
}

double InitialFields::EyInicial( double x )
{
  return 0. ;
}

double InitialFields::EzInicial( double x )
{
  return 0. ;
}

// Define species
Species specie[NSPE] = {
   Species("ions"        ,100.,10,
    0.1,0.001,
    0.,0.5,1.),
   Species("electrons"   ,-1. ,10,
    0.1,0.01 ,
    0.,0.5,1.),
   Species("bg-ions"     ,100.,10,
    0.1,0.001,
    0.5,1.,1.),
   Species("bg-electrons",-1. ,10,
    0.1,0.01 ,
    0.5,1.,1.)
};
