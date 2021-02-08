#include "pic.h"

bool PrintDiagnostics( double t ){
  // Print fields

  std::ofstream file ;

  // Print field values
  for ( int j = 0 ; j < NFIELDS ; j++ ){

    file.open(field_files[j], std::ios::app );
    file << "t = " << t << std::endl ;

    for( int i = 0 ; i < NX ; i++ ){
      file << xx[i] << ' ' << field_var[j][i] << std::endl ;
    }

    file << std::endl ;
    file.close();
  }

  //Print fields
  return true ;
}

void CalculateNewPosVel() {

}

void ComputePIC( const char * dir ){

  // k+1 variables
  std::list<double>* pval1[NSPE][NDIM] ;
  std::list<double>* xval1[NSPE][NDIM] ;

  // Auxiliary vectors for Boris Pusher
  double u[NDIM] , h[NDIM] , s[NDIM];
  double Ek[NDIM] ;

  //Auxiliary variable
  double ql ;
  int num ;

  // Run each time steps of the PIC code
  int i = 0 ;
  for ( double t = 0. ; t < TMAX ; t += dt , i++ ){ // Time steps
    //Step 1 - Compute new positions
    for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species
      num = specie[spe].NumOfPar() ; //Number of Particles of specie.
      ql = specie[spe].ql() ;

      for( int j = 0 ; j < num ; j++ ){ // Particles are in order.
        // Where is the particle located
        int index = int(specie[spe].xval[0][j]/DX) , index1 = index+1;

        //In the case that the particles are in right border
        if ( index == NX ) index1 = 0 ;

        for( int dim = 0 ; dim < NDIM ; dim++ ){ //Scan dimensions

          // Calculate electric field q'E_k
          Ek[dim] = ((xx[index1]-specie[spe].xval[dim][j])*Ex[dim][index] + \
            (specie[spe].xval[dim][j]-xx[index])*Ex[dim][index1])/DX ;
          Ek[dim] *= ql ;

          u[dim] = specie[spe].pval[dim][j] + Ek[dim] ;

          // Calculate h vector (Wiki notation)
          h[dim] = ql*Bx[dim][index] ; //????

          // Calculate s vector (Wiki notation)
          s[dim] = 2*h[dim] ;

          // u = v_{k-1/2}+q'E_k
          u[spe][dim] = specie[spe].xval[dim] ;
        }

        // Finish s calculation
        double aux = 1. ;
        for ( dim = 0 ; dim < NDIM ; dim ++ ) aux += h[dim]*h[dim];
        for ( dim = 0 ; dim < NDIM ; dim ++ ) s[dim] /= aux ;

        
      }

    }

    ComputePosVel();
    if ( i % NDUMP == 0 ){
      if( !PrintDiagnostics(t) ){
        std::cout << "Error in printing the output files." << std::endl ;
        exit(1);
      }
    }
  }

  // Free memory of auxiliary vectors
  for ( int i = 0 ; i < NSPE ; i++ ){
    for ( int j = 0 ; j < NDIM ; j++ ){
      delete[] pval1 ; delete[] xval1 ;

      delete[] u ;
    }
  }

}

void ProduceDiagnostics(){
  std::ofstream file ;
  file.open(output_file,std::ios::app);

  if( !file ){
    std::cout << "ERROR: There was a problem opening the output file." << '\n';
    exit(1) ;
  }

  // Print main parameters
  file << "Number of cores = " << NCOR << std::endl ;
  file << "Number of species = " << NSPE << std::endl ;

  file << "DX = " << dx << std::endl ;
  file << "BOX = " << BOX << std::endl ;

  file << "DT = " << dt << std::endl ;
  file << "TMAX = " << TMAX << std::endl ;

  file.close() ;

  std::cout << "Output file created with success." << '\n';
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
