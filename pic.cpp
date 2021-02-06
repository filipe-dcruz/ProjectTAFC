#include "pic.h"

void ComputePIC(){

  // Run each time steps of the PIC code
  for ( double t = 0. ; t < TMAX ; t += dt ){
    //Step 1
    //Step 4 - Compute new positions
    ComputePosVel();
  }

}

void ProduceDiagnostics(){
  return;

}

void ComputePosVel(){
  return;

}
