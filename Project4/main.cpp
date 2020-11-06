
#include <cstdio>
//#include <omp.h>
#include <cmath>
#include <iostream>

using namespace std;

int * make_spinMatrix(int N){
  int matrix[N*N];

  for (int i = 0; i < N*N; i ++){
    if (rand() % 2 == 1)
      matrix[i] = 1;
    else matrix[i] = -1;
  }

  cout << matrix[1] << matrix[2]<< endl;
  return matrix;
}


int main(int argc, char const *argv[]) {
  int *a;
  a = make_spinMatrix(2);
  cout << *a << endl;
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


  return 0;
}
