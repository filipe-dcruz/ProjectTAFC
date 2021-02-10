/**
    Files responsible by the particle in cell algorith
    @file pic.cpp
    @author Filipe Cruz
**/

#include "pic.h"

/*
  Print the results for a time step in the output files
  @param
    t     : time value of iteration
    pval1 : v_{k+1/2} data
  @return if the print was successful
*/
bool PrintDiagnostics( double t, double* pval1[NSPE][NDIM] ){

  // Offset
  static const double xini[NDIM] = {X0,0.,0.} ;

  static std::ofstream file , file2 ;
  static int numOfParticles , aux , random ;

  std::cout << "DUMP at t = " << t << '\n';

  // Print field values
  for ( int j = 0 ; j < NFIELDS ; j++ ){

    // Open files
    file.open(field_files[j], std::ios::app );
    if ( !file ) return false ;

    // Print time
    file << "t = " << t << std::endl ;

    // Print x and data of files
    for( int i = 0 ; i < NX ; i++ )
      file << xx[i]+xini[0] << ' ' << field_var[j][i] << std::endl ;
    file << std::endl ;

    // close file
    file.close();
  }

  // Print particles values
  for ( int j = 0 ; j < NSPE ; j++ ){

    // Open particle files
    file.open(density_files[j], std::ios::app );
    file2.open(val_files[j], std::ios::app );
    if ( !file || !file2 ) return false ;

    // Print to density files
    file << "t = " << t << std::endl ;
    for( int i = 0 ; i < NX ; i++ )
      file << xx[i]+xini[0] << ' ' << specie[j].density[i] << std::endl ;

    // Print to x and p files
    file2 << "t = " << t << std::endl ;
    numOfParticles = specie[j].NumOfPar() ;
    aux = int(double(numOfParticles)*DUMP_PER) ;

    for ( int k = 0 ; k < numOfParticles ; k++ ){

      // Obtain random to determine which particles to report
      random = rand() % numOfParticles ;

      // Print data to file
      if ( random < aux ){
        for ( int l = 0 ; l < NDIM ; l++ )
          file2 << specie[j].xval[l][k]+xini[l] << ' ' ;
        for ( int l = 0 ; l < NDIM ; l++ )
          file2 << (specie[j].pval[l][k]+pval1[j][l][k])/2. << ' ' ;
        file2 << std::endl ;
      }

    }

    // Close files
    file.close();
    file2.close();
  }

  return true ;
}

/*
  Calculates the phase space data for the next iteration
  @param location to store the new data
*/
void CalculateNewPosVel( double* xval1[NSPE][NDIM], double* pval1[NSPE][NDIM]) {

  // Auxiliary variables
  static double ql ;
  static int num ;

  // Proportions
  static double p1, p2 ;

  // Auxiliary vectors for Boris Pusher
  static double u[NDIM] , h[NDIM] , s[NDIM], ul[NDIM] ;
  static double Ek[NDIM] , auxV[NDIM] ;

  // Step 1 - Compute new positions
  // Scan species
  for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species

    num = specie[spe].NumOfPar() ; //Number of Particles of specie.
    ql = specie[spe].qlvalue() ;

    // Scan number of particles
    for( int j = 0 ; j < num ; j++ ){

      // Where is the particle located
      int index = int(specie[spe].xval[0][j]/dx-0.5) ;
      int index1 = index+1 ;

      //In the case that the particles are in right border
      if ( index1 == NX ) {
        index1 = 0 ;
        p1 = BOX+xx[index1]-specie[spe].xval[0][j] ;
        p2 = specie[spe].xval[0][j]-xx[NX1] ;
      }
      else if ( index == -1 ){
        index = NX1 ;
        p1 = xx[0]-specie[spe].xval[0][j] ;
        p2 = specie[spe].xval[0][j]-xx[NX1]+BOX ;
      }
      else{
        // Calculate electric field q'E_k
        p1 = xx[index1]-specie[spe].xval[0][j] ;
        p2 = specie[spe].xval[0][j]-xx[index] ;
      }

      //Scan dimensions
      for( int dim = 0 ; dim < NDIM ; dim++ ){

        Ek[dim] = (p1*Ex[dim][index] + p2*Ex[dim][index1])*ql ;

        u[dim] = specie[spe].pval[dim][j] + Ek[dim] ;

        // Calculate h vector (Wiki notation)
        h[dim] = ql*(p1*Bx[dim][index]+p2*Bx[dim][index1]) ;

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
        xval1[spe][dim][j] = specie[spe].xval[dim][j]+dt*pval1[spe][dim][j];

      }
    }
  }
}

/*
  Update the phase data with the next iteration data
  @param New data
*/
void UpdateData( double* xval1[NSPE][NDIM], double* pval1[NSPE][NDIM] ){

  static int num ;
  static double ql ;

  // Scan particles species
  for( int spe = 0 ; spe < NSPE ; spe++ ){ // Loop species

    num = specie[spe].NumOfPar() ; //Number of Particles of specie.
    ql = specie[spe].qlvalue() / dx ;

    // In case of borders change x values
    for ( int j = 0 ; j < num ; j++ ){
      if ( xval1[spe][0][j] < 0 )
        xval1[spe][0][j] += BOX ;
      else if ( xval1[spe][0][j] > BOX )
        xval1[spe][0][j] -= BOX ;
    }

    // Scan particles
    for ( int dim = 0 ; dim < NDIM ; dim++ ){
      for ( int j = 1 ; j < num ; j++ ){
        specie[spe].xval[dim][j] = xval1[spe][dim][j] ;
        specie[spe].pval[dim][j] = pval1[spe][dim][j] ;
      }
    }
  }

}

/*
  Calculate a set of constants as function of k once, to increase efficiency
  @param vector to store teh data
*/
void GetFourierVectors( double* res ){

  // Create array for kappa and K
  double kappa , kk ;
  static const double aux3 = dx/(4.*aux1_);

  // Calculate values
  res[0] = 0. ; // consideres for k=0 zero
  for ( int i = 1 ; i < NX ; i++ ){
    kappa = sin(dif1*i) ;
    kk = sin(dif2*i)*sin(dif2*i) ;
    res[i] = aux3*kappa/kk ;
  }

}

/*
  Calculates the next values for the electric field
  @param
    res  : auxiliar value with constants calculated previosly
    Ek   : array to store the field
*/
void FieldSolver( double * res , double * Ek ){

  double auxn , aux2 ;

  // Components of the electric field
  static double Er[NX], Ei[NX] ;
  static double phikr[NX] , phiki[NX] ;
  static double rho[NX] ;

  // Calculate the total density
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

}

/*
  Compute Particle-In-Cell the particle in cell algorith
  @param
    dir  : name of the ouput directory
*/
void ComputePIC( const char * dir ){

  // phase space k+1 variables
  double* pval1[NSPE][NDIM] ;
  double* xval1[NSPE][NDIM] ;

  // E_{k+1}
  double* Ek1 = new double[NX] ;

  // Auxiliar array
  double* res = new double[NX] ;

  //Auxiliary variables
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
  for ( double t = TMIN ; t < TMAX ; t += dt , i++ ){ // Time steps

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

/*
  Function responsable in createing the file with the configuration parameters
*/
void ProduceDiagnostics(){

  // Open output file
  std::ofstream file ;
  file.open(output_file,std::ios::app);

  // Check if file opened
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

  // Close file
  file.close() ;

  std::cout << "Output file created with success." << '\n';
}

/*
  Function responsable in performing the 1D electrostatic PIC code
  @param
    dir   : name of the directory for the output results
*/
void PIC1D( const char * dir = "results/"){

  // Creates files and directory for the results
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

/*
  main function, necessary to start the code. The location of the output
  directory can be given as an argument.
*/
int main(int argc, char const *argv[]) {

  // Run PIC code
  // Checks if a directory name was given and executes PIC simulation
  if ( argc > 1 )
    PIC1D(argv[1]) ;
  else
    PIC1D() ;

	return 0. ;
}
