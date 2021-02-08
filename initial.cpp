#include "initial.h"

double *Ex[NDIM], *Bx[NDIM] ;
double *xx ;

char* name_file ;

char field_files[NFIELDS][MAX_FILE_NAME] = {
  "ex.txt","ey.txt","ez.txt",
  "bx.txt","by.txt","bz.txt"
};

double* field_var[NFIELDS] ;

char output_file[MAX_FILE_NAME] = "output.txt" ;

char density_files[NSPE][MAX_FILE_NAME] ;
char val_files[NSPE][MAX_FILE_NAME] ;

/*
  Initiate parameters for the simulation ;
*/
int CheckParameters( void ){

	// Check invalid values for parameters
  if ( !( NX > 0 ) ){
		std::cout << "--ERROR: INVALID GRID NUMBER" << std::endl ;
		return 1;
	}
  if ( !( BOX > 0. ) ){
		std::cout << "--ERROR: INVALID \"BOX\" VALUE" << std::endl ;
		return 1;
	}
  if ( !( TMAX > TMIN ) ){
		std::cout << "--ERROR: INVALID TIME VALUES" << std::endl ;
		return 1;
	}
  if ( !( CVAL > 0. ) ){
		std::cout << "--ERROR: INVALID COURANT PARAMETER" << std::endl ;
		return 1;
	}
  if ( !( NDUMP > 0. ) ){
		std::cout << "--ERROR: INVALID \"NDUMP\" VALUE" << std::endl ;
		return 1;
	}
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

  //Check parameters
  for( int i = 0 ; i < NSPE ; i++ ){
    if( specie[i].den < 0 ){
      std::cout << "ERROR: Density of specie \"" << specie[i].name
                << "\" is negative." << std::endl ;
      return 1 ;
    }
    if( !specie[i].npar || specie[i].npar > NPARLIMIT ){
      std::cout << "ERROR: Invalid number of particules for specie \""
                << specie[i].name << "\"" << std::endl ;
      return 1 ;
    }
    if( specie[i].xf <= specie[i].x0 ){
      std::cout << "ERROR: Invalid bondaries for specie \"" << specie[i].name
                << "\"" << std::endl ;
      return 1 ;
    }
    if( specie[i].xf > XF || specie[i].x0 < X0  ){
      std::cout << "ERROR: Invalid bondaries for specie \"" << specie[i].name
                << "\"" << std::endl ;
      return 1 ;
    }
    for( int j = 0 ; j < NDIM ; j++ )
      if( specie[i].vth[j] == 0. ){
        std::cout << "ERROR: Thermal velocity is 0. for specie \""
                  << specie[i].name << "\"" << std::endl ;
        return 1 ;
      }

  }
  return 0. ;
}

/*
  Define grids and variables to be used.
*/
void InitialDefinitions( const char * dir ){

  // Seed for random numer
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

  size_t aux = NX * sizeof(double) ;

  // Define grids
  for( int i = 0 ; i < NDIM ; i++ ){
    Ex[i] = (double *) malloc( aux );
    Bx[i] = (double *) malloc( aux );
    field_var[i] = Ex[i] ;
    field_var[i+NDIM] = Bx[i] ;
  }
  xx = (double *) malloc( aux );


  // Declare lists for species
  for ( int i = 0 ; i < NSPE ; i++ ) specie[i].CreateList( NX , dx , dt , X0) ;

  // Create directory for output
  int i = MKDIR , len1 = len + i , len2 = len + add ;

  name_file = new char[len2] ; // Location of directory

  char command[ len1 + add ] = {"mkdir "};
  command[len1] = '/' ;

  for( int j = 0 ; i < len1 ; i++,j++) {
    command[i] = dir[j] ;
    name_file[j] = dir[j];
  }
  name_file[len] = command[len1];

  // Create a new directory for the results
  system(command);
  std::cout << command << '\n';

  // Delete all previous files - case that the directory already existed.
  int len3 = len1 + RM1+RM2 ;
  char command2[ len3 ] = {"rm -r "};
  const char aux2[] = "*.txt" ;

  for ( int i = RM1 , j = 0 ; j < len2 ; i++ , j++ )
    command2[i] = name_file[j] ;
  for ( int i = RM1+len2 , j = 0 ; j < RM2 ; i++ , j++ )
    command2[i] = aux2[j] ;

  // Remove files
  system(command2);
  std::cout << command2 << '\n';

  // Change names of output files to add location
  char * auxN ;
  for (int i = 0 , l = 0 ; i < NFIELDS ; i++ ){
    //Calculate size of file's name
    while(field_files[i][l++]); l-- ;

    // Creates a copy of the name
    auxN = new char[l];
    for( int j = 0 ; j < l; j++ ) auxN[j] = field_files[i][j] ;

    // Update name with location
    for( int j = 0 , k = len2 ; j < l ; j++ , k++ )
      field_files[i][k] = auxN[j] ;
    for( int j = 0 ; j < len2 ; j++ )
      field_files[i][j] = name_file[j];

    delete[] auxN ;
  }

  // For output file of parameters
  int l = 0 ;
  while(output_file[l++]); l-- ;

  auxN = new char[l] ;
  for( int j = 0 ; j < l ; j++ ) auxN[j] = output_file[j] ;

  for( int j = 0 , k = len2 ; j < l ; j++ , k++ )
    output_file[k] = auxN[j] ;
  for( int j = 0 ; j < len2 ; j++ )
    output_file[j] = name_file[j] ;

  delete[] auxN ;

  // Output files for Particles
  for( int i = 0 ; i < NSPE ; i++ ){
    int m ;
    for ( int j = 0 ; j < len2 ; j++ ){
      density_files[i][j] = name_file[j] ;
      val_files[i][j] = name_file[j] ;
    }

    l = 0 ;
    while(specie[i].name[l++]); l-- ;
    m = len2 + l ;

    for ( int j = len2 , k = 0 ; k < l ; j++ , k++ ){
      density_files[i][j] = specie[i].name[k] ;
      val_files[i][j] = specie[i].name[k] ;
    }

    l = 0 ; while(density_name[l++]); l-- ;
    for( int j = m, k = 0 ; k < l ; j++ , k++ )
      density_files[i][j] = density_name[k] ;

    l = 0 ; while(val_name[l++]); l-- ;
    for( int j = m, k = 0 ; k < l ; j++ , k++ )
      val_files[i][j] = val_name[k] ;

  }

}

/*
  Define parameters with the setup of the simulation ;
*/
void DefineInitialValues( void ){

  static double x = dx/2. ; // X0 removed for convinience
  static itr it1, it2 ;
  int j , k ;

  // initiate grids for electromagnetic fields
  for( int i = 0 ; i < NX ; i++ , x += dx ){
    Ex[0][i] = InitialFields::ExInicial( x ) ;
    Ex[1][i] = InitialFields::EyInicial( x ) ;
    Ex[2][i] = InitialFields::EzInicial( x ) ;
    Bx[0][i] = InitialFields::BxInicial( x ) ;
    Bx[1][i] = InitialFields::ByInicial( x ) ;
    Bx[2][i] = InitialFields::BzInicial( x ) ;
    xx[i] = x ;
  }

  // Species
  for ( int i = 0 ; i < NSPE ; i++ ){

    // Auxiliary - distances between particles
    double auxDx  = ( specie[i].xf - specie[i].x0 ) / specie[i].NumOfPar() ;
    double auxDx2 = auxDx / 2. ;

    x = specie[i].x0+auxDx2 ;

    // Random generator for each specie for the Maxwellian distribution
    for ( j = 0 ; j < NDIM ; j++ ){
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(
        specie[i].v0[j],specie[i].vth[j]);

      it1 = specie[i].xval[j]->begin() ;
      it2 = specie[i].pval[j]->begin() ;

      // Uniform distribution of particles
      for( ; it1 != specie[i].xval[j]->end() ; it1++ , it2++ , x += auxDx ){
        *it1=x ;                      // Add position
        *it2=distribution(generator); // Add velocity
      }

      // Quick way to deal with other components
      x = 0. ;
      auxDx = 0. ;
    }

    // Calculate density
    double qcval = specie[i].qcvalue() / dx ;
    it1 = specie[i].xval[0]->begin() ;
    it2 = specie[i].xval[0]->end() ;

    double aux ;

    for ( j = 0 ; j < NX ; j++ ) specie[i].density[j] = 0. ;

    // Left particles
    while( (aux = xx[1]-(*it1)) > 0. && it1 != it2 ){
      specie[i].density[0] += dx-(*it1-xx[0]) ;
      specie[i].density[1] += dx-aux ;
      it1++;
    }

    // Middle particles
    for ( j = 1 , k = 2 ; j < NX1 ; j++ , k++ ){
      while( (aux = xx[k]-(*it1)) > 0. && it1 != it2 ){
        specie[i].density[j] += dx-(*it1-xx[j]);
        specie[i].density[k] += dx-aux ;
        it1++;
      }
    //  std::cin.get() ;
    }

    // Right particles
    while( (aux = BOX-(*it1)) > 0. && it1 != it2 ){ // Right side particles
      specie[i].density[NX1] += dx-(*it1-xx[NX1]) ;
      specie[i].density[0] += dx-aux ;
      it1++;
    }

    //Multiple by cloud charge
    for ( j = 0 ; j < NX ; j++ ) specie[i].density[j] *= qcval ;

  }
}

void FinalDeclarations( void ){
  // free name of directory
  free(name_file);

  // Free memory of grids
  for ( int i = 0 ; i < NDIM ; i++ ){
    free(Ex[i]);
    free(Bx[i]);
  }
  free(xx);
}
