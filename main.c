/**
    @file pic1d.c
    @author Filipe Cruz
		@descrition Code that performes particle in cell code in 1D
*/

#include <stdlib.h>
#include <stdio.h>

#include "input.h" //Input configuration
#include "initial.h" //Input configuration
#include "pic.h" //PIC code

int main(int argc, char const *argv[]) {

	// initiate parameters
	printf("\n--STARTING PROGRAM--\n\n");
	InitialDeclarations();

	// initiate remaining paramentes with the configuration
	printf("Initiate variables...\n");
	DefineInitialValues();

	// Compute results
	printf("Compute results...\n");

	// Print paramenters of simulation
	printf("Produce output...\n");
	FinalDeclarations() ;

	// Finish program
	printf("\n--FINISHING PROGRAM--\n\n");

	return 0. ;
}
