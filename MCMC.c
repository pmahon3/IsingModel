#include "MCMC.h"


/* randomizeLattice either randomly initilizes the lattice or assigns all spins up
   p is the pointer to the lattice
   n is the number of rows
   m is the number of columns
   nonrandom, if 1 assigns all spins up, randomly assigns spins otherwise
*/
void randomizeLattice(arrPtr p, int n, int m, int nonrandom) {

  if ( nonrandom ){
    for ( int i = 0 ; i < n ; i++ ){
      for ( int j = 0 ; j < m ; j++ ){
	*(( p + i * m ) + j ) = 1;
      }
    }
    return;
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {

      double r = (double) genrand();
      if ( r < 0.5 ) *(( p + i * m) + j) = 1;
      else *(( p + i*m ) + j ) = -1;
    }
  }
}

/* totalHamiltonian computes the total Hamiltonian energy of the system
   p is the pointer to the lattice
   n is the number of rows in the lattice
   m is the number of columns in the lattice
   J is the coupling constant
*/
int totalHamiltonian(arrPtr p, int n, int m, int J){
  int total = 0;
  int spinij, spinijright, spinijbottom;
  for ( int i = 0 ; i < n ; i++ ){
    for ( int j = 0 ; j < m ; j++ ){
      spinij = *( p + i*m + j);
      spinijright = *( p + i*m + modulo( j + 1, m ) );
      spinijbottom = *( p + modulo( i * m + j + m, n*m ));
      total += spinij * ( spinijright + spinijbottom);
    }
  }
  return (total * -1 * J);
}

/* totalMagnetization computes the total magentization of the system
   p is the pointer to the lattice
   n is the number of rows in the lattice
   m is the number of columns in the lattice
*/
int totalMagnetization(arrPtr p, int n, int m){
  int total = 0;
  for ( int i = 0 ; i < n ; i++ ){
    for ( int j = 0 ; j < m; j++ ){
      total +=  *(p + i * m + j);
    }
  }
  return total;
}

/* transitionProbability specifices the transition function and computes the probability of changing the state
   dE is Hamiltonian energy of the system given the proposed change
   t is the temperature of the system
*/
double transitionProbability(int dE, double t) {
  double exponent =  (double) (-1 * dE) / (double) (KB*t) ;
  double numer = powf(E, exponent);  
  double denom = 1 + powf(E, exponent);
  double result = (double) numer / (double) denom;
  return result;
}

/* dE computes the change in Hamiltonian energy if the given spin is flipped.
   p is the pointer to the lattice
   i is the row of the spin to be flipped
   j is the column of the spin to be flipped
   n is the number of rows
   m is the number of columns
   J is the coupling constant
 */
int dESpin(arrPtr p, int i, int j, int n, int m, int J) {

  
  int topshift, bottomshift, rightshift, leftshift, top, bottom, right, left, sum, spinij;

  spinij = *( p + i * m  + j);

  topshift = modulo( (i * m + j - m), ( n*m ));
  bottomshift = modulo( (i * m + j + m), ( n*m ));
  rightshift = ( i * m ) + modulo( j + 1 , m );
  leftshift = ( i * m ) + modulo( j - 1, m );

  top = *(p + topshift);
  bottom = *(p + bottomshift);
  right = *(p + rightshift);
  left = *(p + leftshift);
  
  return -1*J * -2 * spinij * (top + bottom + right + left );
}

/* flipSpin flips the given spin
   p is the pointer to the lattice
   i is the row of the spin to be flipped
   j is the column of the spin to be flipped
   n is the number of rows
   m is the number of columns
 */
void flipSpin( arrPtr p, int i, int j, int n, int m){
  int *spinij = ( p + i * m + j);
  *spinij = -1 * *(spinij);
}

/* printLattice prints the lattice to the terminal
   p is the pointer to the lattice
   n is the number of rows
   m is the number of columns
 */
void printLattice(arrPtr p, int n, int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      printf("%3d", *( p + i * m  + j ));
    }
    printf("\n");
  }
  printf("\n");
}

/* modulo implenets the modulus function for integer inputs mapping to the natural numbers
   a is the number modulated
   n is the modulus
 */
int modulo(int a, int n){
  int mod = a % n;
  if ( a < 0 ){
    mod += n;
  }
  return mod;  
}
