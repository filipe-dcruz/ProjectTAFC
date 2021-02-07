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

void PIC1D(){

  // initiate parameters
  std::cout << "\n--STARTING PROGRAM--\n" << std::endl ;
  InitialDefinitions();

  // initiate remaining paramentes with the configuration
  std::cout << "Initiating variables..." << std::endl ;
  DefineInitialValues();

	// Compute results
	std::cout << "Done\nStarting simulation..." << std::endl ;
	ComputePIC() ;

	// Print paramenters of simulation
	std::cout << "Simulation done.\nProducing diagnostics..." << std::endl ;
	ProduceDiagnostics();

	// Finish program
	std::cout << "Diagnostics done\n\n--FINISHING PROGRAM--\n" << std::endl ;
	FinalDeclarations() ;

  return;
}

int main(int argc, char const *argv[]) {

	//Run PIC code
	PIC1D() ;

	return 0. ;

}
