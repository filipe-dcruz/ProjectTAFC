/**
    @file pic1d.c
    @author Filipe Cruz
		@descrition Code that performes particle in cell code in 1D plasmas
*/

#include <stdio.h>

#include "input.h"    // Input configuration
#include "initial.h"  // Input configuration
#include "pic.h"      // PIC code

int main(int argc, char const *argv[]) {

	// initiate parameters
	printf("\n--STARTING PROGRAM--\n\n");
	InitialDeclarations();

	// initiate remaining paramentes with the configuration
	printf("Initiating variables...");
	DefineInitialValues();

	// Compute results
	printf("Done\nStarting simulation...");
	ComputePIC() ;

	// Print paramenters of simulation
	printf("Simulation done.\nProducing diagnostics...\n");
	ProduceDiagnostics();

	// Finish program
	printf("Diagnostics done\n\n--FINISHING PROGRAM--\n\n");
	FinalDeclarations() ;

	return 0. ;
}
