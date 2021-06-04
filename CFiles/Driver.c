#include "Graph.h"
#include "MT19937.h"
#include "AuxillaryFunctions.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#define NTHREADS 6

void *mainfun(void *x){

  // Get thread id
  int tid = *((int *) x);

  // Create a results file for the thread

  char buf[0x100];
  snprintf(buf, sizeof(buf), "Thread%d.csv", tid); 

  
  FILE* results  = fopen(buf, "w");

  sgenrand(time(0));

  int tsteps = 2000;
  int burnin = 100;
  int size = 3000;

  int nreplicates = 3;

  // Mean node degree
  int m0 = 2;

  // Graph initialization
  graph* g = createGraph(size, m0, 0, 0);

  // Run simulations

  double beta = 0.05;
  double delta = 0.9;
  double score = beta / delta * g->lambda;

  while( delta >= 0.01){

    g->phi = beta;
    g->eta = delta;

    for ( int i = 0; i < nreplicates ; i++ ){
      printf("Thread %d\nLambda: %f\nBeta: %f\nDelta: %f\nScore: %f\nSimulation %d/%d\n\n", tid, g->lambda, g->phi, g->eta, score, i+1, nreplicates);
      simEpidemic(g, burnin, tsteps);

      printf("Thread %d\nAverage ratio of infected: %f\nSimulation %d/%d\n\n", tid, g->avginfected, i+1, nreplicates );
      fprintf(results, "%f, %f,\n", score, g->avginfected );
    }
    delta -= 0.01;
    score = beta / delta * g->lambda;
  }

  fclose(results);

  printf("Thread %d complete.\n", tid);

  return NULL;
}

int main(){

  // Parallelization code sourced from gribblelab.org/CBootCamp/A2_Parallel_Programming_in_C.html

  clock_t start, end;
  double cpu_time_used;
  start = clock();
  
  pthread_t threads[NTHREADS];
  int thread_args[NTHREADS];
  int rc, i;

  /* spawn the threads */
  for (i=0; i<NTHREADS; ++i)
    {
      thread_args[i] = i;
      printf("spawning thread %d\n", i);
      rc = pthread_create(&threads[i], NULL, mainfun, (void *) &thread_args[i]);
    }

  /* wait for threads to finish */
  for (i=0; i<NTHREADS; ++i) {
    rc = pthread_join(threads[i], NULL);
  }

  end = clock();

  cpu_time_used = ((double) (end - start )) / CLOCKS_PER_SEC;
  printf("\n\nTOTAL EXCECUTION TIME:\t%f", cpu_time_used);
  return 0;
}
