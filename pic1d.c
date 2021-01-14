// 1-d PIC code to solve plasma two-stream instability problem.

#include <stdlib.h>
#include <stdio.h>

const char INPUT_DEFAULT_NAME[] = "input.txt" ;
const char INPUT_DEFAULT_TYPE[] = ".txt" ;

int main(int argc, char const *argv[]) {

	// Initial constants
	const char *input_name ; //Name of file

  printf("\n-----START OF PROGRAM-----\n");

	// Read input data
	if( argc < 2 ){ //No args
		input_name = INPUT_DEFAULT_NAME ;
	}
	else if ( argc == 2 ){
		//Check file used
		if ( argv[1] == ".txt")
			input_name = INPUT_DEFAULT_NAME ;
		else {
			printf("ERROR: Invalid name for input file\n");
			return 1 ;
		}
	}
	else {
		printf("ERROR: Invalid number of arguments\n");
		return 1 ;
	}

	printf("-----END OF PROGRAM-----\n\n");

	return 0;
}
