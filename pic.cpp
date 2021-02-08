#include "pic.h"

bool PrintDiagnostics( double t ){

  std::cout << "DUMP at t = " << t << '\n';

  static std::ofstream file , file2 ;

  static int aux , aux2 , random ;

  itr it1[NDIM] , it2[NDIM];

  // Print field values
  for ( int j = 0 ; j < NFIELDS ; j++ ){

    file.open(field_files[j], std::ios::app );
    if ( !file ) return false ;

    file << "t = " << t << std::endl ;

    for( int i = 0 ; i < NX ; i++ )
      file << xx[i] << ' ' << field_var[j][i] << std::endl ;

    file << std::endl ;
    file.close();
  }

  // Print particles values
  for ( int j = 0 ; j < NSPE ; j++ ){

    file.open(density_files[j], std::ios::app );
    file2.open(val_files[j], std::ios::app );

    if ( !file || !file2 ) return false ;

    // Print to density files
    file << "t = " << t << std::endl ;
    for( int i = 0 ; i < NX ; i++ )
      file << xx[i] << ' ' << specie[j].density[i] << std::endl ;

    // Print to x and p files
    file2 << "t = " << t << std::endl ;
    aux = specie[j].NumOfPar() ;
    aux2 = int(double(aux)*DUMP_PER) ;

    for ( int k = 0 ; k < NDIM ; k++ ){
      it1[k] = specie[j].xval[k]->begin();
      it2[k] = specie[j].pval[k]->begin();
    }

    while( it1[0] != specie[j].xval[0]->end() ){

      // Obtain random to determine which particles to report
      random = rand() % aux ;

      if ( random < aux2 ){
        for ( int k = 0 ; k < NDIM ; k++ ) file2 << *it1[k] << ' ' ;
        for ( int k = 0 ; k < NDIM ; k++ ) file2 << *it2[k] << ' ' ;
        file2 << std::endl ;
      }

      // increment iterators
      for ( int k = 0 ; k < NDIM ; k++ ) {
        it1[k]++;
        it2[k]++;
      }

    }

    file.close();
    file2.close();
  }

  return true ;
}

void CalculateNewPosVel() {

}

void ComputePIC( const char * dir ){

  // k+1 variables
  double* pval1[NSPE][NDIM] ;
  double* xval1[NSPE][NDIM] ;

  itr it1[NDIM] , it2[NDIM] , it11[NDIM] , it21[NDIM] ;

  // Auxiliary vectors for Boris Pusher
  double u[NDIM] , h[NDIM] , s[NDIM], ul[NDIM] ;
  double Ek[NDIM] , auxV1[NDIM] , auxV2[NDIM] ;

  //Auxiliary variable
  double ql ;
  int num ;

  // Proportions
  double p1, p2 ;

  // Allocate memory
  for( int i = 0 ; i < NSPE ; i++ )
    for( int j = 0 ; j < NDIM ; j++ ){
      pval1[i][j] = new double [specie[i].NumOfPar()];
      xval1[i][j] = new double [specie[i].NumOfPar()];
    }


  // Run each time steps of the PIC code
  int i = 0 ;
  for ( double t = 0. ; t < TMAX ; t += dt , i++ ){ // Time steps

    std::cout << "t = " << t << '\n';

    // Print Diagnostics
    if ( i % NDUMP == 0 ){
      if( !PrintDiagnostics(t) ){
        std::cout << "Error in printing the output files." << std::endl ;
        exit(1);
      }
    }

    //Step 1 - Compute new positions
    for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species
      num = specie[spe].NumOfPar() ; //Number of Particles of specie.
      ql = specie[spe].qlvalue()/dx ;

      std::cout << "Updated positions" << '\n';
      for ( int dim = 0 ; dim < NDIM ; dim++ ){
        it1[dim] = specie[spe].xval[dim]->begin() ;
        it2[dim] = specie[spe].pval[dim]->begin() ;
      }

      std::cout << "j" << '\n';
      for( int j = 0 ; j < num ; j++ ){ // Particles are in order.
        std::cout << "k" << '\n';

        // Where is the particle located
        int index = int(*it1[0]/dx) ;
        int index1 = index+1 ;

        //In the case that the particles are in right border
        if ( index == NX ) index1 = 0 ;

        for( int dim = 0 ; dim < NDIM ; dim++ ){ //Scan dimensions
          std::cout << index << '\n';
          // Calculate electric field q'E_k
          p1 = xx[index1]-(*it1[dim]) ;
          p2 = (*it1[dim])-xx[index] ;
          Ek[dim] = (p1*Ex[dim][index] + p2*Ex[dim][index1]) ;
          Ek[dim] *= ql ;

          u[dim] = *it2[dim] + Ek[dim] ;

          // Calculate h vector (Wiki notation)
          h[dim] = ql*(p1*Bx[dim][index]+p2*Bx[dim][index1]) ; //????

          // Calculate s vector (Wiki notation)
          s[dim] = 2*h[dim] ;
        }
        std::cout << "Updated"<<j<<' '<<spe << '\n';

        // Finish s calculation
        double aux = 1. ;
        for ( int dim = 0 ; dim < NDIM ; dim ++ ) aux += h[dim]*h[dim];
        for ( int dim = 0 ; dim < NDIM ; dim ++ ) ul[dim] /= aux ;

        // Calculate u'
        CrossProduct(u,h,auxV1);
        for ( int dim = 0 ; dim < NDIM ; dim ++ ) auxV1[dim] += u[dim] ;

        CrossProduct(auxV1,s,auxV2);
        for ( int dim = 0 ; dim < NDIM ; dim ++ ) {
          ul[dim] += u[dim] + auxV2[dim];

          // Calculate v_{l+1/2}
          pval1[spe][dim][j] = ul[dim] + Ek[dim];

          //Calculate x_{k+1}
          xval1[spe][dim][j] = (*it1[dim])+dt*pval1[spe][dim][j];

          // Incremente interators
          it1[dim]++;
          it2[dim]++;
        }
      }
    }

    for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species

      std::cout << "Updated positions" << '\n';

      // Set iterators to Upgrade data
      for ( int dim = 0 ; dim < NDIM ; dim++ ){
        it11[dim] = specie[spe].xval[dim]->begin() ;
        it21[dim] = specie[spe].pval[dim]->begin() ;
        it1[dim] = it11[dim]++ ;
        it2[dim] = it21[dim]++ ;

        // update the first elements j = 0
        (*it1[dim]) = xval1[spe][dim][0] ;
        (*it2[dim]) = pval1[spe][dim][0] ;
      }

      std::cout << "Updated int" << '\n';

      // Upgrade data
      for ( int j = 1 , k = 0 ; j < num ; j++ , k++ )
      {
        // Update values
        std::cout << spe << ' ' << j<< ' '<< num << ' ' << xval1[spe][0][j] << '\n';
        for ( int dim = 0 ; dim < NDIM ; dim++ ){
          (*it11[dim]) = xval1[spe][dim][j] ;
          (*it21[dim]) = pval1[spe][dim][j] ;
        }
        std::cout << "oi"<< j <<' '<< k << '\n';

        // Swap necessary particles
        if ( xval1[spe][0][j] < xval1[spe][0][k] ){
          std::cout << "o" << '\n';
          for ( int dim = 0 ; dim < NDIM ; dim++ ){
            std::swap(*it1[dim],*(it11[dim])) ;
            std::swap(*it2[dim],*(it21[dim])) ;
            it1[dim]++;it2[dim]++;
          }
          j++; k++;
        }
        std::cout << "oi2" << '\n';

        //Increment iterators
        if ( it1[0] != specie[spe].xval[0]->end() ){ // If it there is no more elements
          for ( int dim = 0 ; dim < NDIM ; dim++ ){
            it1[dim]++ ;
            it2[dim]++ ;
          }
        }
      }

    }
  }

  // Free memory
  for( int i = 0 ; i < NSPE ; i++ )
    for( int j = 0 ; j < NDIM ; j++ ){
      delete[] pval1[i][j] ;
      delete[] xval1[i][j] ;
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
