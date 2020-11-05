#include <cstdio>
#include <omp.h>
#include <cmath>
#include <iostream>

using namespace std;

int main(int argc, char const *argv[]) {
  /*
  int n_spins, **spin_matrix, mcs;
  long idum;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
  double Delta_E[5];
  double Beta =
  // Calculates all exp(Beta * Delta_E)
  for (int i = 0; i < 5; i ++)
      Delta_E[i] = exp(Beta * 4*(i-2));
  */

  for (int i = 0; i < 10; i++)
      cout << rand() % 2 << endl; //return int between 0 and 1
  return 0;
}
