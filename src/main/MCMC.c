#include "MCMC.h"

//--2D-ISING-MODEL-FUNCTIONS------------------------------------------------------------------------------------------//

/*  Description:
        randomizeLattice either randomly initializes the lattice or assigns all spins up

    Inputs:
        p   The pointer to the lattice
        n   The number of rows
        m   The number of columns
        nonrandom   If 1 assigns all spins up, randomly assigns spins otherwise

    Output:
        Modifies array p in place
*/
void randomize_lattice(arrPtr p, int n, int m, int nonrandom) {
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

/*  Description:
        totalHamiltonian computes the total Hamiltonian energy of the system.

    Inputs:
        p           The pointer to the lattice
        n           The number of rows
        m           The number of columns
        J           The isotropic coupling constant

    Returns:
        The Hamiltonian of the the lattice.
*/
int total_hamiltonian(arrPtr p, int n, int m, int J){
    int total = 0;
    int spin_ij, spin_ij_right, spin_ij_bottom;
    for ( int i = 0 ; i < n ; i++ ){
        for ( int j = 0 ; j < m ; j++ ){
            spin_ij = *( p + i*m + j);
            spin_ij_right = *( p + i*m + modulo( j + 1, m ) );
            spin_ij_bottom = *( p + modulo( i * m + j + m, n*m ));
            total += spin_ij * ( spin_ij_right + spin_ij_bottom);
        }
    }
    return (total * -1 * J);
}

/*  Description:
        totalMagnetization computes the total Hamiltonian energy of the system.

    Inputs:
        p           The pointer to the lattice
        n           The number of rows
        m           The number of columns

    Returns:
        The total magnetization for the lattice.
*/
int total_magnetization(arrPtr p, int n, int m){
    int total = 0;
    for ( int i = 0 ; i < n ; i++ ){
        for ( int j = 0 ; j < m; j++ ){
            total +=  *(p + i * m + j);
        }
    }
  return total;
}

/*  Description:
        dEspin computes the change in Hamiltonian energy if the given spin is flipped.

    Inputs:
        p   A pointer to the lattice
        i   The row of the spin to be flipped
        j   The column of the spin to be flipped
        n   The number of rows
        m   The number of columns
        J   The isotropic coupling constant

    Returns:
        The change in energy for the proposed spin flip
 */
int dE_spin(arrPtr p, int i, int j, int n, int m, int J) {
    int top_shift, bottom_shift, right_shift, left_shift, top, bottom, right, left, sum, spin_ij;

    spin_ij = *( p + i * m  + j);

    top_shift = modulo( (i * m + j - m), ( n*m ));
    bottom_shift = modulo( (i * m + j + m), ( n*m ));
    right_shift = ( i * m ) + modulo( j + 1 , m );
    left_shift = ( i * m ) + modulo( j - 1, m );

    top = *(p + top_shift);
    bottom = *(p + bottom_shift);
    right = *(p + right_shift);
    left = *(p + left_shift);

    return -1*J * -2 * spin_ij * (top + bottom + right + left );
}

/*  Description:
        flipSpin flips the given spin

    Inputs:
        p   A pointer to the lattice
        i   The row of the spin to be flipped
        j   The column of the spin to be flipped
        n   The number of rows
        m   The number of columns

    Output:
        Modifies array p in place.
 */
void flip_spin( arrPtr p, int i, int j, int n, int m){
    int *spin_ij = ( p + i * m + j);
    *spin_ij = -1 * *(spin_ij);
}

/*  Description:
        printLattice prints the lattice to the terminal.

    Inputs:
        p   A pointer to the lattice
        n   The number of rows
        m   The number of columns

    Output:
        Prints to standard out.

 */
void print_lattice(arrPtr p, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%3d", *( p + i * m  + j ));
        }
        printf("\n");
    }
    printf("\n");
}

//--MCMC-FUNCTIONS----------------------------------------------------------------------------------------------------//

/*  Description:

*/

void mcmcStep(arrPtr p, int n, int m){
    return;
}


/*  Description:
        transitionProbability specifies the transition function for the system and computes the probability of a state
        change

    Inputs:
        dE  Hamiltonian energy differential for the proposed state transition
        t   The temperature of the system

    Returns:
        The probability of a state change occurring.
*/

double transition_probability(int dE, double t) {
    double exponent =  (double) (-1 * dE) / (double) (KB*t);
    double numer = powf(E, exponent);
    double denom = 1 + powf(E, exponent);
    double result = (double) numer / (double) denom;
    return result;
}

/*  Description:
        Modulo implements the modulus function for integer inputs mapping to the natural numbers

    Inputs:
        a   The number to be modulated
        n   The modulus

    Outputs:
        a mod n
 */
int modulo(int a, int n){
    int mod = a % n;
    if ( a < 0 ){
        mod += n;
    }
    return mod;
}
