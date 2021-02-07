#include "initial.h"

double *Ex, *Ey, *Ez ;
double *Bx, *By, *Bz ;
double *xx ;

char* name_file ;

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
  // Check species names
  for( int i = 0 ; i < NSPE ; i++ )
    for( int j = i+1 , k ; j < NSPE ; j++ ){
      for ( int k = 0 ; k < NAME_LIMIT ; k++)
        if( specie[i].name[k] != specie[j].name[k] ) break ;
      if( k == NAME_LIMIT ){
        std::cout << "--ERROR: Species " << i << " and " << j
                  << " have the same name." << std::endl ;
    		return 1;
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
    if( specie[i].vth == 0. ){
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

  int len = 0 ;
  while ( dir[len++] );

  int add = ( dir[len-2] == '/' ? 0 : 1 ) ;

  // Check parameters
  if ( CheckParameters() ){
    std::cout << "Program will terminate..." << std::endl ;
    exit(1) ;
  }

  size_t aux = NX * sizeof(double) ;

  // Define grids
  Ex = (double *) malloc( aux );
  Ey = (double *) malloc( aux );
  Ez = (double *) malloc( aux );
  Bx = (double *) malloc( aux );
  By = (double *) malloc( aux );
  Bz = (double *) malloc( aux );
  xx = (double *) malloc( aux );

  // Species
  for ( int i = 0 ; i < NSPE ; i++ ) specie[i].CreateList( NX , dx ) ;

  // Create directory for output
  int i = MKDIR , len1 = len + i , len2 = len + add ;

  name_file = new char[len2] ;

  char command[ len1 + add ] = {"mkdir "};
  command[len1] = '/' ;

  for( int j = 0 ; i < len1 ; i++,j++) {
    command[i] = dir[j] ;
    name_file[j] = dir[j];
  }
  name_file[len1] = dir[len1];

  system(command);

  // Delete all previous files - case that the directory already existed.
  int len3 = len1 + RM1+RM2 ;
  char command2[ len3 ] = {"rm "};
  const char aux[] = "*.txt" ;
  
  for ( int i = 3 , j = 0 ; j < len2 ; i++ , j++ ) command2[i] = name_file[j] ;
  for ( int i = len2 , j = 0 ; j < 5 ; i++ , j++ ) command2[i] = aux[j] ;

  system(command2);

}

/*
  Define parameters with the setup of the simulation ;
*/
void DefineInitialValues(void ){

  double x = dx/2. ;

  // initiate grids
  for( int i = 0 ; i < NX ; i++ , x += dx ){
    Ex[i] = InitialFields::ExInicial( x ) ;
    Ey[i] = InitialFields::EyInicial( x ) ;
    Ez[i] = InitialFields::EzInicial( x ) ;
    Bx[i] = InitialFields::BxInicial( x ) ;
    By[i] = InitialFields::ByInicial( x ) ;
    Bz[i] = InitialFields::BzInicial( x ) ;
    xx[i] = x ;
  }

  // Species
  for ( int i = 0 ; i < NSPE ; i++ ){

    // Auxiliary
    double auxDx  = ( specie[i].xf - specie[i].x0 ) / specie[i].NumOfPar() ;
    double auxDx2 = auxDx / 2. ;

    x = specie[i].x0+auxDx2 ;

    // Random generator for each specie
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(specie[i].v0,specie[i].vth);

    //Ignore initial cells without particles
    //int ini = int((x-X0)/dx) ;
    //for( int j = 0 ; j < ini ; j++ ){
    //  specie[i].xpos[j] = NULL ;
    //  specie[i].ppos[j] = NULL ;
    //}

    itr it1 = specie[i].xval->begin() ;
    itr it2 = specie[i].pval->begin() ;

    // Uniform distribution of particles
    for( ; it1 != specie[i].xval->end() ; it1++ , it2++ , x += auxDx ){
      *it1=x ;
      *it2=distribution(generator);
      //if( xx[ini] < x ) {
      //  specie[i].xpos[ini]   = it1 ;
      //  specie[i].ppos[ini++] = it2 ;
      //}
    }

    // Ignore final cells without particules
    //for( int j = ini ; j < NX ; j++ ){
    //  specie[i].xpos[j] = NULL ;
    //  specie[i].ppos[j] = NULL ;
    //}
  }
}

void FinalDeclarations( void ){
  // Free memory
  free(Ex); free(Ey); free(Ez);
  free(Bx); free(By); free(Bz);
  free(xx);

  // free name of directory
  free(name_file);
}
