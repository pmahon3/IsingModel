#include "MCMC.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

int main() {

    // Simulation parameters and values
    int n_chains, n_burn_in, n_iters, n_thin, n, Mag, Ham, J, n_temps, mag_sum, ham_sum, temp, i, j, dE;

    double *tempPtr, spacing, current_temp, magnetism_ensemble_avg, hamiltonian_ensemble_avg, prob, rand;
    
    // Output folder for final results
    char out_folder[50] = "../resources/out/";
    struct stat st = {0};
    if (stat(out_folder, &st) == -1) {
        mkdir(out_folder, 0700);
    }

    // Results csv
    FILE *resultsfile = fopen(strcat(out_folder,"results_file.csv"), "w");;
    
    // Seed random generator and set values
    sgenrand(time(0));
    
    n_burn_in = 3000;
    n_iters = 3000;
    n_thin = 10;
    n = 64;
    region_size = 4
    J = 1;
    n_temps = 100;
    
    // Declare and populate array for simulation temperatures
    tempPtr = calloc(n_temps, sizeof(double ) );
    
    spacing = (double)5/(double)n_temps;
    
    for( int i = 0 ; i < n_temps ; i++ ){
        double *pos = tempPtr + i;
        *pos = spacing * (double)(i+1);
    }


    // Perform MCMC simulation for each temperature
    for( int temp = 0 ; temp < n_temps ; temp++ ){

        // Declare array for lattice
        int arr_ptr [n][n];
        double current_temp = *(tempPtr + temp);

        // Output files for results of burn in period and sampling period
        char burn_file_name[60];
        char sample_file_name[60];

        sprintf(burn_file_name, "%sburn_file_%f.csv", out_folder, current_temp);
        sprintf(sample_file_name, "%ssample_file_%f.csv", out_folder, current_temp);

        FILE *burn_file = fopen(burn_file_name, "w");
        FILE *sample_file = fopen(sample_file_name, "w");

        printf("Simulating temperature %f ...\n", current_temp);

        // Randomize initial configuration of lattice
        randomize_lattice( (int *) arr_ptr, n, n, 0);

        // Print initial lattice
        printf( "Initial lattice...\n\n");
        printLattice( (int *) arr_ptr, n, n);

        Mag = total_magnetization( (int *) arr_ptr, n, n );
        Ham = total_hamiltonian( (int *) arr_ptr, n, n, J );

        // Burn in phase
        printf("Burn in...\n");

        // Compute n_burn_in MCMC steps
        for ( int burn = 0 ; burn < n_burn_in ; burn++ ){
            // Explore n-by-n dimensional space
            for ( int dim = 0 ; dim < n * n ; dim++ ){
                // Randomly select one dimension and calculate change in Hamiltonian if state flipped
                i = floor(genrand() * (double) n);
                j = floor(genrand() * (double) n);
                dE = dE_spin( (int *) arr_ptr, i, j, n, n, J );

                // Compute transition probability and update state
                prob = transition_probability(dE, current_temp);
                if ( prob > genrand() ){
                  flip_spin( (int *) arr_ptr, i, j, n, n);
                  Mag = Mag + 2 * *((int *) arr_ptr + i * n + j);
                  Ham = Ham + dE;
                }
            }

            // Accept every n_thin-th MCMC step
            if ( modulo( burn, n_thin) == 0){
                // Save actual magnetism and hamiltonian
                fprintf(burn_file, "%d,%d,%d,\n", burn/n_thin, Mag, Ham);

                // Compute regional magnetism and hamiltonian
                if (region_size > 1){
                    n_regions = ( n * n ) / ( region_size * region_size );
                    regions_per_row = n / region_size;

                    for ( int i = 0; i < n_regions; i++ ){
                        row_ptr = (arr_ptr + ((i * region_size) / n) * region_size) * n;
                        col_ptr = (i % regions_per_row) * region_size;
                        region_ptr = row_ptr + col_ptr
                    }
                }
            }
        }
    }

    fclose(burn_file);
    printf("Burn in finished...\n\n");

    // Print lattice after burn in period
    //printLattice((int *) arr_ptr, n, n);

    // Sampling period
    printf("Sampling...\n");

    // Initialize statistical measures
    mag_sum = 0;
    ham_sum = 0;
    hamiltonian_ensemble_avg = 0;
    magnetism_ensemble_avg = 0;

    // Compute n_iters MCMC steps
    for ( int step = 0 ; step < n_iters ; step++ ){

      // Explore n-by-n dimensional space
      for( int dim = 0 ; dim < n * n ; dim++ ){
        // Randomly select one dimension and calculate change in Hamiltonian if state flipped
        i = floor(genrand() * (double) n);
        j = floor(genrand() * (double) n);
        dE = dE_spin( (int *) arr_ptr, i, j, n, n, J);

        // Compute transition probability and update state
        prob = transition_probability(dE, current_temp);
        if ( prob > genrand() ){
          flip_spin( (int *) arr_ptr, i, j, n, n);
          Mag = Mag + 2 * *((int *) arr_ptr + i * n + j);
          Ham = Ham + dE ;
	    }
      }

      // Accept n_thin-th accept MCMC step
      if( modulo(step, n_thin) == 0 ){
        fprintf(sample_file, "%d,%d,%d,\n", step/n_thin, Mag, Ham);
        ham_sum += Ham;
        mag_sum += abs(Mag);
      }
    }

    // Compute summary statistics
    hamiltonian_ensemble_avg =  (double) ham_sum / (double) ( n_iters / n_thin );
    magnetism_ensemble_avg = (double) mag_sum / (double) ( n_iters / n_thin );

    // Write out to results file
    fprintf(resultsfile, "%f,%f,%f,%f,%f,\n", current_temp, hamiltonian_ensemble_avg, magnetism_ensemble_avg, 
    hamiltonian_ensemble_avg / (double) (n*n), magnetism_ensemble_avg / (double) (n*n));

    // Print final lattice
    print_lattice( (int *) arr_ptr, n, n);
    
    printf("Complete\n\n");
    fclose(sample_file);  
  }
  fclose(resultsfile);
  return 0;
}

// Initialize

// Burn in

// Run

// ---- //

// Transition Probability

// MCMC step



