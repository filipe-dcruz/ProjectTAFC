#include "pic.h"

bool PrintDiagnostics( double t, double* pval1[NSPE][NDIM] ){

  static const double xini[NDIM] = {X0,0.,0.} ;

  std::cout << "DUMP at t = " << t << '\n';

  static std::ofstream file , file2 ;

  static int aux , aux2 , random ;

  itr it1[NDIM] , it2[NDIM] ;
  double* it3[NDIM] ;

  // Print field values
  for ( int j = 0 ; j < NFIELDS ; j++ ){

    file.open(field_files[j], std::ios::app );
    if ( !file ) return false ;

    file << "t = " << t << std::endl ;

    for( int i = 0 ; i < NX ; i++ )
      file << xx[i]-xini[0] << ' ' << field_var[j][i] << std::endl ;

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
      file << xx[i]-xini[0] << ' ' << specie[j].density[i] << std::endl ;

    // Print to x and p files
    file2 << "t = " << t << std::endl ;
    aux = specie[j].NumOfPar() ;
    aux2 = int(double(aux)*DUMP_PER) ;

    for ( int k = 0 ; k < NDIM ; k++ ){
      it1[k] = specie[j].xval[k]->begin();
      it2[k] = specie[j].pval[k]->begin();
      it3[k] = pval1[j][k];
    }

    while( it1[0] != specie[j].xval[0]->end() ){

      // Obtain random to determine which particles to report
      random = rand() % aux ;

      if ( random < aux2 ){
        for ( int k = 0 ; k < NDIM ; k++ )
          file2 << (*it1[k])-xini[k] << ' ' ;
        for ( int k = 0 ; k < NDIM ; k++ )
          file2 << ((*it2[k])+(*it3[k]))/2. << ' ' ;
        file2 << std::endl ;
      }

      // increment iterators
      for ( int k = 0 ; k < NDIM ; k++ ) {
        it1[k]++;
        it2[k]++;
        it3[k]++;
      }

    }

    file.close();
    file2.close();
  }

  return true ;
}

void CalculateNewPosVel( double* xval1[NSPE][NDIM], double* pval1[NSPE][NDIM]) {

  //Auxiliary variable
  static double ql ;
  static int num ;

  // Proportions
  static double p1, p2 ;

  static itr it1[NDIM] , it2[NDIM] , it11[NDIM] , it21[NDIM] ;

  // Auxiliary vectors for Boris Pusher
  static double u[NDIM] , h[NDIM] , s[NDIM], ul[NDIM] ;
  static double Ek[NDIM] , auxV[NDIM] ;

  //Step 1 - Compute new positions
  for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species
    num = specie[spe].NumOfPar() ; //Number of Particles of specie.
    ql = specie[spe].qlvalue() / dx ;

    for ( int dim = 0 ; dim < NDIM ; dim++ ){
      it1[dim] = specie[spe].xval[dim]->begin() ;
      it2[dim] = specie[spe].pval[dim]->begin() ;
    }

    for( int j = 0 ; j < num ; j++ ){ // Particles are in order.

      // Where is the particle located
      int index = int(*it1[0]/dx) ;
      if (index < 0) std::cout << "Here" << j << ' ' << *it1[0] <<' ' << num << ' '<<spe<< '\n';
      int index1 = index+1 ;

      //In the case that the particles are in right border
      if ( index == NX ) index1 = 0 ;

      for( int dim = 0 ; dim < NDIM ; dim++ ){ //Scan dimensions
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

      // Finish s calculation
      double aux = 1. ;
      for ( int dim = 0 ; dim < NDIM ; dim ++ ) aux += h[dim]*h[dim];
      for ( int dim = 0 ; dim < NDIM ; dim ++ ) s[dim] /= aux ;

      // Calculate u'
      CrossProduct(u,h,auxV);
      for ( int dim = 0 ; dim < NDIM ; dim ++ ) auxV[dim] += u[dim] ;

      CrossProduct(auxV,s,ul);
      for ( int dim = 0 ; dim < NDIM ; dim ++ ) {
        ul[dim] += u[dim] ;

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

}

void UpdateData( double* xval1[NSPE][NDIM], double* pval1[NSPE][NDIM] ){

  static const double offset[NDIM] = {BOX,0.,0.};

  static int num ;
  static double ql ;

  static itr it1[NDIM] , it2[NDIM] , it11[NDIM] , it21[NDIM] ;

  for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species
    num = specie[spe].NumOfPar() ; //Number of Particles of specie.
    ql = specie[spe].qlvalue() / dx ;

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

    // Upgrade data
    for ( int j = 1 , k = 0 ; j < num ; j++ , k++ ){
      // Update values
      for ( int dim = 0 ; dim < NDIM ; dim++ ){
        (*it11[dim]) = xval1[spe][dim][j] ;
        (*it21[dim]) = pval1[spe][dim][j] ;
      }

      // Swap necessary particles
      if ( xval1[spe][0][j] < xval1[spe][0][k] ){
        for ( int dim = 0 ; dim < NDIM ; dim++ ){
          std::swap(*it1[dim],*it11[dim]) ;
          std::swap(*it2[dim],*it21[dim]) ;
          it1[dim]++; it11[dim]++;
          it2[dim]++; it21[dim]++;
        }
        j++; k++;
      }

      //Increment iterators
      if ( j < num ){ // If there is no more elements
        for ( int dim = 0 ; dim < NDIM ; dim++ ){
          it1[dim]++ ;
          it2[dim]++ ;
          it21[dim]++;
          it11[dim]++;
        }
      }

    }

    // Auxiliary
    for ( int dim = 0 ; dim < NDIM ; dim++ ){
      it1[dim] = specie[spe].xval[dim]->begin() ;
      it2[dim] = specie[spe].pval[dim]->begin() ;
      it11[dim] = specie[spe].xval[dim]->end() ;
      it21[dim] = specie[spe].pval[dim]->end() ;
      it11[dim]-- ;
      it21[dim]-- ;
    }

    bool t1 = (*(it1[0]) < 0.) , t2 = (*(it11[0]) > BOX);

    // Swap border particles
    // If first passed to last
    if( t1 && t2 ){ // Swap both ends
      for ( int dim = 0 ; dim < NDIM ; dim++ ){
        *(it1[dim]) += offset[dim] ;
        *(it11[dim]) -= offset[dim] ;
        std::swap(*(it1[dim]),*(it11[dim])) ;
        std::swap(*(it2[dim]),*(it21[dim])) ;
      }
    }
    else if ( t1 ){
      for ( int dim = 0 ; dim < NDIM ; dim++ ){
        *(it1[dim]) += offset[dim] ;
        specie[spe].xval[dim]->push_back(*(it1[dim])) ; // Add Last
        specie[spe].pval[dim]->push_back(*(it2[dim])) ; // Add Last
        specie[spe].xval[dim]->pop_front() ;            // Remove first
        specie[spe].pval[dim]->pop_front() ;            // Remove first
      }
    }
    else if ( t2 ){
      for ( int dim = 0 ; dim < NDIM ; dim++ ){
        *(it11[dim]) -= offset[dim] ;
        specie[spe].xval[dim]->push_front(*(it11[dim])) ; // Add First
        specie[spe].pval[dim]->push_front(*(it21[dim])) ; // Add First
        specie[spe].xval[dim]->pop_back() ;            // Remove Last
        specie[spe].pval[dim]->pop_back() ;            // Remove Last
      }
    }

  }

}

// Get vectors that will be used to calculate the field
void GetFourierVectors( double* res ){

  // Create array for kappa and K
  double kappa , kk ;

  double aux3 = dx/(4.*aux1_);

  // Calculate values
  res[0] = 0. ;
  for ( int i = 1 ; i < NX ; i++ ){
    kappa = sin(dif1*i) ;
    kk = sin(dif2*i)*sin(dif2*i) ;
    res[i] = aux3*kappa/kk ;
  }

}

void FieldSolver( double * res , double * Ek ){

  double auxn , aux2 ;

  // Components of the electric field
  static double Er[NX], Ei[NX] ;
  static double phikr[NX] , phiki[NX] ;
  static double rho[NX] ;

  // Calculate total density
  for ( int i = 0 ; i < NX ; i++ )
    rho[i] = specie[0].density[i] ;
  for ( int j = 1 ; j < NSPE ; j++ )
    for( int i = 0 ; i < NX ; i++ )
      rho[i] += specie[j].density[i] ;

  // Calculate rho_k
  for ( int n = 0 ; n < NX ; n++ ){
    phikr[n] = 0. ;
    phiki[n] = 0. ;
    auxn = aux_*n ;

    for ( int i = 0 ; i < NX ; i++ ){ // get each k
      aux2 = auxn*xx[i];
      phikr[n] += rho[i]*cos(aux2) ;
      phiki[n] -= rho[i]*sin(aux2) ;
    }

    // Calculate poterntial in k
    phikr[n] *= res[n] ; phiki[n] *= res[n] ;

    // Multiply by -i
    Er[n] = phiki[n] ;
    Ei[n] = -phikr[n] ;
  }

  // IFFT - real part
  for ( int i = 0 ; i < NX ; i++ ){
    Ek[i] = 0. ;
    auxn = aux_*xx[i] ;
    for ( int n = 0 ; n < NX ; n++ ){
      aux2 = auxn*n ;
      Ek[i] += (Er[n]*cos(aux2)-Ei[n]*sin(aux2)) ;
    }
    Ek[i] /= BOX ;
  }

  /*Ek1[0] = _dx2*(phi[NX1]-phi[1]) ;
  Ek1[NX1] = _dx2*(phi[NX2]-phi[0]) ;
  for ( int i = 1 , j = 2 , k = 0 ; i < NX1 ; i++ , k++, l++ )
    for Ek1[i] = _dx2*(phi[k]-phi[j]) ;
  */
}

/*
  Compute Particle-In-Cell
*/
void ComputePIC( const char * dir ){

  // k+1 variables
  double* pval1[NSPE][NDIM] ;
  double* xval1[NSPE][NDIM] ;

  // E_{k+1}
  double* res = new double[NX] ;
  double* Ek1 = new double[NX] ;

  //Auxiliary variable
  double ql ;
  int num ;

  // Calculate kappa and kk vectors for the field solver
  GetFourierVectors( res ) ;

  // Allocate memory for k+1 and k+1/2 data
  for( int i = 0 ; i < NSPE ; i++ )
    for( int j = 0 ; j < NDIM ; j++ ){
      pval1[i][j] = new double [specie[i].NumOfPar()];
      xval1[i][j] = new double [specie[i].NumOfPar()];
    }

  // Run each time steps of the PIC code
  int i = 0 ;
  for ( double t = 0. ; t < TMAX ; t += dt , i++ ){ // Time steps

    // Boris pusher
    CalculateNewPosVel(xval1,pval1) ;

    // Print Diagnostics
    // Print after particles to do average of velocity
    if ( i % NDUMP == 0 ){
      if( !PrintDiagnostics(t,pval1) ){
        std::cout << "Error in printing the output files." << std::endl ;
        exit(1);
      }
    }

    // Electric field solver
    FieldSolver(res,Ek1);

    // Update data
    UpdateData(xval1,pval1) ;

    // Update density from the updated positions
    CalculateTheDensity() ;

    // Update electric field
    for ( int j = 0 ; j < NX ; j++ ) Ex[0][j] = Ek1[j] ;
  }

  // Free memory
  for( i = 0 ; i < NSPE ; i++ )
    for( int j = 0 ; j < NDIM ; j++ ){
      delete[] pval1[i][j] ;
      delete[] xval1[i][j] ;
    }

  delete[] Ek1 ;
  delete[] res ;

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

void PIC1D( const char * dir = "results/"){

  std::cout << "dir" << dir << '\n';

  // initiate parameters
  std::cout << "\n--STARTING PROGRAM--\n" << std::endl ;
  InitialDefinitions( dir );

  // initiate remaining paramentes with the configuration
  std::cout << "Initiating variables..." << std::endl ;
  DefineInitialValues();

	// Compute results
	std::cout << "Done\nStarting simulation..." << std::endl ;
	ProduceDiagnostics();
	ComputePIC(dir) ;

	// Print paramenters of simulation
	std::cout << "Simulation done." << std::endl ;

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
