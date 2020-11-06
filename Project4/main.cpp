#include <cstdio>
#include <omp.h>
#include <cmath>
#include <iostream>

using namespace std;

int*[] make_matrix(matrix, N){
  for (int i = 0; i < N*N; i ++){
    if (rand() % 2 == 1)
      matrix[i] = 1;
    else matrix[i] = -1;
  }
}


int main(int argc, char const *argv[]) {

  int* matrix[] = make_matrix(N);
  int *matrix2 = make_matrix(N_2);


  /*
  int n_spins, *spin_matrix, mcs;
  long idum;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
  double Delta_E[5];
  double Beta =
  // Calculates all exp(Beta * Delta_E)
  for (int i = 0; i < 5; i ++)
      Delta_E[i] = exp(Beta * 4*(i-2));
  */

  initialize_matrix(spin_matrix);

  return 0;
}
