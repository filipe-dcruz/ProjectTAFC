#include "pic.h"

bool PrintDiagnostics( double t ){
  // Print fields

  std::ofstream file ;

  // Print field values
  for (int j = 0 ; j < NFIELDS ; j++ ){

    file.open(field_files[j], std::ios::app );
    file << "t = " << t << std::endl ;

    for( int i = 0 ; i < NX ; i++ ){
      file << xx[i] << ' ' << field_var[j] << std::endl ;
    }

    file << std::endl ;
    file.close();
  }

  //Print fields
  return true ;
}

void ComputePIC( const char * dir ){

  // Run each time steps of the PIC code
  int i = 0 ;
  for ( double t = 0. ; t < TMAX ; t += dt , i++ ){
    //Step 1
    //Step 4 - Compute new positions
    ComputePosVel();
    if (i % NDUMP){
      if( !PrintDiagnostics(t) ){
        std::cout << "Error in printing the output files." << std::endl ;
        exit(1);
      }
    }
  }

}

void ProduceDiagnostics(){
  return;
}

void ComputePosVel(){
  return;

}

void PIC1D( const char * dir = "results/"){

  // initiate parameters
  std::cout << "\n--STARTING PROGRAM--\n" << std::endl ;
  InitialDefinitions( dir );

  // initiate remaining paramentes with the configuration
  std::cout << "Initiating variables..." << std::endl ;
  DefineInitialValues();

	// Compute results
	std::cout << "Done\nStarting simulation..." << std::endl ;
	ComputePIC(dir) ;

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
  if ( argc > 1 )
    PIC1D(argv[1]) ;
  else
    PIC1D() ;

	return 0. ;

}
