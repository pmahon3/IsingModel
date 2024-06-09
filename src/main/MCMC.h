#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* MCMC parameters */
#define E 2.718281828
#define KB 1

typedef int nchains;
typedef int nburnin;
typedef int niters;
typedef int nthin;
typedef int nrows;
typedef int ncols;
typedef int J;
typedef int ntemps;
typedef int *arrPtr;

/* Mersenne twister parameters */
/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

/* Random generator seed */
typedef unsigned long seed;


/* MCMC functions */
void randomize_lattice(arrPtr p, int n, int m, int nonrandom);
int total_hamiltonian(arrPtr p, int n, int m, int J);
int total_magnetization(arrPtr p, int n, int m);
double transition_probability(int dE, double t);
void flip_spin(arrPtr p, int i, int j, int n, int m);
int dE_spin(arrPtr p, int i, int j, int n, int m, int J);
int modulo( int a, int n);
void print_lattice(arrPtr p, int n, int m);

/* Mersenne twister functions */ 
void sgenrand(seed s);
double genrand();
