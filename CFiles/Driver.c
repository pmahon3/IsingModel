#include "MCMC.h"
#include <unistd.h>

int main() {

  // Simulation parameters and values
  int nchains, nburnin, niters, nthin, nrows, ncols, Mag, Ham, J, ntemps, magsum, hamsum, temp, i, j, dE;

  double *tempPtr, spacing, currenttemp, magensemble, hamensemble, prob, rand;

  // Output for final results
  FILE *resultsfile = fopen("ResultsFile.csv", "w");;

  // Seed random generator and set values
  sgenrand(time(0));

  nburnin = 3000;
  niters = 3000;
  nthin = 10;
  nrows = 30;
  ncols = 30;
  J = 1;
  ntemps = 100;

  // Declare and populate array for simulation temperatures
  tempPtr = calloc(ntemps, sizeof(double ) );

  spacing = (double)5/(double)ntemps;
  
  for( int n = 0 ; n < ntemps ; n++ ){
    double *pos = tempPtr + n;
    *pos = spacing * (double)(n+1);
  }


  // Perform MCMC simulation for each temperature
  for( int temp = 0 ; temp < ntemps ; temp++ ){

    // Declare array for lattice
    int arrPtr [nrows][ncols];
    currenttemp = *(tempPtr + temp);

    // Output files for results of burn in period and sampling period
    char burnbuf[0x100], samplebuf[0x100];
    snprintf(burnbuf, sizeof(burnbuf), "BurnFile%f.csv", currenttemp);
    snprintf(samplebuf, sizeof(samplebuf), "SampleFile%f.csv", currenttemp);
    FILE *burnfile = fopen(burnbuf, "w");
    FILE *samplefile = fopen(samplebuf, "w");
    
    printf("Simulating temperature %f ...\n", currenttemp);

    // Randomize initial configuration of lattice
    randomizeLattice( (int *) arrPtr, nrows, ncols, 0);

    // Print initial lattice
    printf( "Initial lattice...\n\n");
    printLattice( (int *) arrPtr, nrows, ncols);

    Mag = totalMagnetization( (int *) arrPtr, nrows, ncols );
    Ham = totalHamiltonian( (int *) arrPtr, nrows, ncols, J );

    // Burn in phase
    printf("Burn in...\n");

    // Compute nburnin MCMC steps
    for( int burn = 0 ; burn < nburnin ; burn++ ){


      // Explore nrows-by-ncols dimensional space
      for ( int dim = 0 ; dim < nrows * ncols ; dim++ ){

	// Randomly select one dimension and calculate change in Hamiltonian if state flipped
	i = floor(genrand() * (double) nrows);
	j = floor(genrand() * (double) ncols);
	dE = dESpin( (int *) arrPtr, i, j, nrows, ncols, J );

	// Compute transition probability and update state
	prob = transitionProbability(dE, currenttemp);
	if ( prob > genrand() ){
	  flipSpin( (int *) arrPtr, i, j, nrows, ncols);
	  Mag = Mag + 2 * *((int *) arrPtr + i * ncols + j);
	  Ham = Ham + dE;	
	}
      }

      // Accept every nthin-th MCMC step
      if ( modulo( burn, nthin) == 0) fprintf(burnfile, "%d,%d,%d,\n", burn/nthin, Mag, Ham);
      
    }

    fclose(burnfile);
    printf("Burn in finished...\n\n");

    // Print lattice after burn in period
    //printLattice((int *) arrPtr, nrows, ncols);

    // Sampling period
    printf("Sampling...\n");

    // Initialize statistical measures
    magsum = 0;
    hamsum = 0;
    hamensemble = 0;
    magensemble = 0;

    // Compute niters MCMC steps
    for ( int step = 0 ; step < niters ; step++ ){

      // Explore nrows-by-ncols dimensional space
      for( int dim = 0 ; dim < nrows * ncols ; dim++ ){
	// Randomly select one dimension and calculate change in Hamiltonian if state flipped
	i = floor(genrand() * (double) nrows);
	j = floor(genrand() * (double) ncols);
	dE = dESpin( (int *) arrPtr, i, j, nrows, ncols, J);

	// Compute transition probability and update state
	prob = transitionProbability(dE, currenttemp);
	if ( prob > genrand() ){
	  flipSpin( (int *) arrPtr, i, j, nrows, ncols);
	  Mag = Mag + 2 * *((int *) arrPtr + i * ncols + j);
	  Ham = Ham + dE ;
	}
      }

      // Accept nthin-th accept MCMC step
      if( modulo(step, nthin) == 0 ){
	fprintf(samplefile, "%d,%d,%d,\n", step/nthin, Mag, Ham);
	hamsum += Ham;
	magsum += abs(Mag);
      }
    }

    // Compute summary statistics
    hamensemble =  (double) hamsum / (double) ( niters / nthin );
    magensemble = (double) magsum / (double) ( niters / nthin );

    // Write out to results file
    fprintf(resultsfile, "%f,%f,%f,%f,%f,\n", currenttemp, hamensemble, magensemble, hamensemble / (double) (nrows*ncols), magensemble / (double) (nrows*ncols));

    // Print final lattice
    printLattice( (int *) arrPtr, nrows, ncols);
    
    printf("Complete\n\n");
    fclose(samplefile);  
  }
  fclose(resultsfile);
  return 0;
}
