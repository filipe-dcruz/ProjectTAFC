/**
    Declaration and definition of multiple calculated parameters and functions
    responsable for the initiation of the variables and objects.
    @file initial.cpp
    @author Filipe Cruz
**/

#include "initial.h"

// Declaration of previous undeclared paramentes in header file
double *Ex[NDIM], *Bx[NDIM] ;
double *xx  ;
char* name_file ;

// Names for the electric field files
char field_files[NFIELDS][MAX_FILE_NAME] = {
  "ex.txt","ey.txt","ez.txt",
  "bx.txt","by.txt","bz.txt"
};

double* field_var[NFIELDS] ;

char output_file[MAX_FILE_NAME] = "output.txt" ;

char density_files[NSPE][MAX_FILE_NAME] ;
char val_files[NSPE][MAX_FILE_NAME] ;

/*
  Check if the parameter of the simulation are valid quantities.
  @return 1 if error and zero in case of success
*/
int CheckParameters( void ){

	// Check invalid values for parameters
  if ( !( NX > 0 ) ){ // If number of cells is negative
		std::cout << "--ERROR: INVALID GRID NUMBER" << std::endl ;
		return 1;
	}
  if ( !( BOX > 0. ) ){ // If box size is negative
		std::cout << "--ERROR: INVALID \"BOX\" VALUE" << std::endl ;
		return 1;
	}
  if ( !( TMAX > TMIN ) ){ // If swapped time limits
		std::cout << "--ERROR: INVALID TIME VALUES" << std::endl ;
		return 1;
	}
  if ( !( CVAL > 0. ) ){ // If negative Corant number
		std::cout << "--ERROR: INVALID COURANT PARAMETER" << std::endl ;
		return 1;
	}
  if ( !( NDUMP > 0. ) ){ // If negative dumping frequency
		std::cout << "--ERROR: INVALID \"NDUMP\" VALUE" << std::endl ;
		return 1;
	}
  // If Courant number is bigger than 1, it may lead to resolution problems
	if ( CVAL > 1. )
		std::cout << "--WARNING: Courant number is higher than 1.0" << std::endl ;

  // Check species

  // Check if species names are repeated
  for( int i = 0 ; i < NSPE ; i++ ){
    for( int j = i+1 , k ; j < NSPE ; j++ ){

      for ( k = 0 ; k < NAME_LIMIT && specie[i].name[k] && specie[j].name[k] ; k++ )
        if( specie[i].name[k] != specie[j].name[k] ) break ;

      if( k == NAME_LIMIT || (!specie[i].name[k] && !specie[j].name[k]) ){
        std::cout << "--ERROR: Species " << i << " and " << j
                  << " have the same name." << std::endl ;
    		return 1;
      }
    }
  }

  //Check parametersof Species
  for( int i = 0 ; i < NSPE ; i++ ){
    // If negative density
    if( specie[i].den < 0 ){
      std::cout << "ERROR: Density of specie \"" << specie[i].name
                << "\" is negative." << std::endl ;
      return 1 ;
    }
    // If invalid number of particles per cell
    if( specie[i].npar <= 0 || specie[i].npar > NPARLIMIT ){
      std::cout << "ERROR: Invalid number of particules for specie \""
                << specie[i].name << "\"" << std::endl ;
      return 1 ;
    }
    // If swapped boundaries
    if( specie[i].xf <= specie[i].x0 ){
      std::cout << "ERROR: Invalid bondaries for specie \"" << specie[i].name
                << "\"" << std::endl ;
      return 1 ;
    }
    // If boundaries are out of the box
    if( specie[i].xf > XF || specie[i].x0 < X0  ){
      std::cout << "ERROR: Invalid bondaries for specie \"" << specie[i].name
                << "\"" << std::endl ;
      return 1 ;
    }
    // If invalid thermal velocities
    for( int j = 0 ; j < NDIM ; j++ )
      if( specie[i].vth[j] == 0. ){
        std::cout << "ERROR: Thermal velocity is 0. for specie \""
                  << specie[i].name << "\"" << std::endl ;
        return 1 ;
      }

  }

  // Success
  return 0. ;

}

/*
  Define grids and variables to be used and creates names for output files
  @ param
    dir - name of folder to store the output results
*/
void InitialDefinitions( const char * dir ){

  // Seed for random number generator
  srand (time(NULL));

  // Get length of directory name
  int len = 0 ;
  while ( dir[len++] );
  len-- ;

  // See if '/' was included
  int add = ( dir[len-1] == '/' ? 0 : 1 ) ;

  // Check if parameters values are okay
  if ( CheckParameters() ){
    std::cout << "Program will terminate..." << std::endl ;
    exit(1) ;
  }

  // Necessary memory for the grids
  size_t aux = NX * sizeof(double) ;

  // Allocate electromagnetic grids
  for( int i = 0 , i1 = NDIM; i < NDIM ; i++ , i1++ ){
    Ex[i] = (double *) malloc( aux );
    Bx[i] = (double *) malloc( aux );
    field_var[i] = Ex[i] ;
    field_var[i1] = Bx[i] ;
  }

  // Allocate cells position
  xx = (double *) malloc( aux );

  // Declare lists for each specie
  for ( int i = 0 ; i < NSPE ; i++ )
    specie[i].CreateList( NX , dx , dt , X0) ;

  // Create directory for output files
  // Auxiliary lengths for the directory and commands
  int i = MKDIR , len1 = len + i , len2 = len + add ;

  name_file = new char[len2] ; // Location of directory

  // Command necessary
  char command[ len1 + add ] = {"mkdir "};
  command[len1] = '/' ;

  // Add directory to command
  for( int j = 0 ; i < len1 ; i++,j++) {
    command[i] = dir[j] ;
    name_file[j] = dir[j];
  }
  name_file[len] = command[len1];

  // Create a new directory for the results
  system(command);
  std::cout << command << '\n';

  // Delete all previous files - case that the directory already existed.
  // Auxiliar lengths
  int len3 = len1 + RM1+RM2 ;

  // Command rm
  char command2[ len3 ] = {"rm -r "};
  const char aux2[] = "*.txt" ;

  // Add directory to command
  for ( int i = RM1 , j = 0 ; j < len2 ; i++ , j++ )
    command2[i] = name_file[j] ;
  for ( int i = RM1+len2 , j = 0 ; j < RM2 ; i++ , j++ )
    command2[i] = aux2[j] ;

  // Remove files
  system(command2);
  std::cout << command2 << '\n';

  // Change names of output files to add location of directory
  char * auxN ;
  for (int i = 0 , l = 0 ; i < NFIELDS ; i++ ){

    //Calculate size of file's name
    while(field_files[i][l++]); l-- ;

    // Creates a copy of the name
    auxN = new char[l];
    for( int j = 0 ; j < l; j++ ) auxN[j] = field_files[i][j] ;

    // Update name with location of directory
    for( int j = 0 , k = len2 ; j < l ; j++ , k++ )
      field_files[i][k] = auxN[j] ;
    for( int j = 0 ; j < len2 ; j++ )
      field_files[i][j] = name_file[j];

    // Free temporary copy
    delete[] auxN ;

  }

  // For output files of the parameters of the simulation
  // Calculate length of name
  int l = 0 ;
  while(output_file[l++]); l-- ;

  // Create copy of name
  auxN = new char[l] ;
  for( int j = 0 ; j < l ; j++ ) auxN[j] = output_file[j] ;

  // Add directory location
  for( int j = 0 , k = len2 ; j < l ; j++ , k++ )
    output_file[k] = auxN[j] ;
  for( int j = 0 ; j < len2 ; j++ )
    output_file[j] = name_file[j] ;

  delete[] auxN ;

  // Output files for Particles
  for( int i = 0 ; i < NSPE ; i++ ){
    int m ;

    // Add location of directory
    for ( int j = 0 ; j < len2 ; j++ ){
      density_files[i][j] = name_file[j] ;
      val_files[i][j] = name_file[j] ;
    }

    // Calculate length of species name
    l = 0 ;
    while(specie[i].name[l++]); l-- ;
    m = len2 + l ;

    // Change the output files with the species name
    for ( int j = len2 , k = 0 ; k < l ; j++ , k++ ){
      density_files[i][j] = specie[i].name[k] ;
      val_files[i][j] = specie[i].name[k] ;
    }

    // Change for the density files
    l = 0 ; while(density_name[l++]); l-- ;
    for( int j = m, k = 0 ; k < l ; j++ , k++ )
      density_files[i][j] = density_name[k] ;

    // Change for the phase space files
    l = 0 ; while(val_name[l++]); l-- ;
    for( int j = m, k = 0 ; k < l ; j++ , k++ )
      val_files[i][j] = val_name[k] ;

  }

}

/*
  Auxiliar function to calculate the density from the
  positions of the particles. Uses the general approach.
*/
void CalculateTheDensity( void ){
  // Particles are not sorted, so we cannot use other efficient ways to
  // calculate the density

  //Auxiliar variables
  static double qcval , aux ;
  int index , index1 ;

  // Scan through the number of species
  for ( int i = 0 ; i < NSPE ; i++ ){

    // Total number of particles
    int numOfParticles = specie[i].NumOfPar() ;

    // Calculate density
    qcval = specie[i].qcvalue() ;

    // initiate densities to zero
    for ( int j = 0 ; j < NX ; j++ ) specie[i].density[j] = 0. ;

    // Scan the particles
    for ( int k = 0 ; k < numOfParticles ; k++ ){
      // Get index of particle
      index = int(specie[i].xval[0][k]/dx-0.5) ; // I need to look at the center
      index1 = index+1 ;

      if( index1 == NX ) { // Particles at the right
        specie[i].density[NX1] += BOX+xx[0]-specie[i].xval[0][k] ;
        specie[i].density[0] += specie[i].xval[0][k]-xx[NX1] ;
      }
      else if ( index == -1 ){ // Particles at the left
        specie[i].density[NX1] += xx[0]-specie[i].xval[0][k] ;
        specie[i].density[0] += specie[i].xval[0][k]-xx[NX1]+BOX;
      }
      else{
        // Add contribution of cloud (CIC model)
        specie[i].density[index] += xx[index1]-specie[i].xval[0][k] ;
        specie[i].density[index1] += specie[i].xval[0][k]-xx[index] ;
      }

    }

    //Multiply by cloud charge
    for ( int j = 0 ; j < NX ; j++ ) specie[i].density[j] *= qcval ;
  }

}

/*
  Define parameters with the setup of the simulation
*/
void DefineInitialValues( void ){

  static double x = dx/2. ; // X0 removed for convinience

  // Initiate grids for electromagnetic fields and cell's positions
  for( int i = 0 ; i < NX ; i++ , x += dx ){
    Ex[0][i] = InitialFields::ExInicial( x ) ;
    Ex[1][i] = InitialFields::EyInicial( x ) ;
    Ex[2][i] = InitialFields::EzInicial( x ) ;
    Bx[0][i] = InitialFields::BxInicial( x ) ;
    Bx[1][i] = InitialFields::ByInicial( x ) ;
    Bx[2][i] = InitialFields::BzInicial( x ) ;
    xx[i] = x ;
  }

  // Scan through the species
  for ( int i = 0 ; i < NSPE ; i++ ){

    // Total number of particles
    int numOfParticles = specie[i].NumOfPar() ;

    // Auxiliary - distances between particles
    double auxDx  = ( specie[i].xf - specie[i].x0 ) / specie[i].NumOfPar() ;
    double auxDx2 = auxDx / 2. ;

    // Position of particle
    x = specie[i].x0+auxDx2 ;

    // Random generator for each orientation for the Maxwellian distribution
    for ( int j = 0 ; j < NDIM ; j++ ){

      // Start gaussian generator
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(
        specie[i].v0[j],specie[i].vth[j]);

      // Uniform distribution of particles positions
      for( int k = 0 ; k < numOfParticles ; k++, x += auxDx ){
        specie[i].xval[j][k] = x ;                       // Add position
        specie[i].pval[j][k] = distribution(generator) ; // Add velocity
      }

      // Quick way to deal with other components
      x = 0. ;
      auxDx = 0. ;
    }
  }

  // Calculate the density from the particle positions
  CalculateTheDensity() ;

}

/**
  Release temporary allocated memory of the grids
**/
void FinalDeclarations( void ){

  // Free name of directory
  free(name_file);

  // Free memory of grids
  for ( int i = 0 ; i < NDIM ; i++ ){
    free(Ex[i]);
    free(Bx[i]);
  }

  // Free positions of the cells
  free(xx);
}
