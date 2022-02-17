#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef MCMC{

	
	int nIterations;
	int nBurnIn;
	int nThin;
	int nChains;
	
	int T; 
	int size;

	struct Lattice* currentLattice;
	struct Lattice* previousLattice;


	void (*initialize)(const struct MCMC*, struct Lattice*, int nIterations, int size, int nBurnIn, int nThin, int nChains, int T);
	void (*runChains)(const struct MCMC*);
	void (*plotChains)(const struct MCMC*);
	bool (*transitionProbability)(int dE);
};

void initialize(const struct MCMC* sim, struct Lattice* initialLattice, int nIterations, int size, int nBurnIn, int nThins, int nChains, int T) {

	sim.currentLattice = *initialLattice;
	sim.previousLattice = *initialLattice;
	sim.nIterations = nIterations;
	sim.size = size;
	sim.nBurnIn = nBurnIn;
	sim.nThins = nThins;
	sim.nChains = nChains;
	sim.T = T;
};

void runChains(const struct MCMC* sim) {

	for (int i = 0 ; i < sim.nBurnIn ; i++) {
		for (int j = 0; j < sim.nSize; j++) {

			

			int row = rand() % sim.nRows;
			int col = rand() % sim.nCols; 

			int dE = sim.currentLattice.dESpin(row, col);

			if (sim.transitionProbability(sim.T, dE)) {

				sim.previousLattice = *sim.currentLattice;
				sim.currentLattice.flipSpin(row, col);
			}
		}
	}
}

bool transitionProbability( int T, int dE) {
	
	srand((unsigned) time(&t));
		
	double prob = exp((double)(-1 * dE / T)) / (1 + exp(-1 * dE / T)) );

	double num = (double)rand() / (double)RAND_MAX;

	if (prob < num) return true;
	else return false;
}