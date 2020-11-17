
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "isingModel.hpp"

using namespace std;


int main(int argc, char const *argv[]) {
  //double Ti = stod(argv[2]);
  //double Tf = stod(argv[3]);
  //double dT = stod(argv[4]);



  int L = atoi(argv[1]);

  double Ti = 2.0;
  double Tf = 2.35;
  double dT = 0.00125;
  int T_length = 32;
  double T_array [T_length];

  for (int i = 0; i < T_length; i++){
    T_array[i] = Ti + dT*i;
    cout << T_array[i] << endl;
  }

  ofstream Lfile;
  Lfile.open("Observables_" + to_string(L) + ".csv");
  Lfile <<  "T, <E>, <M>, Cv, chi" << endl;

  double start;
  double end;
  // 283 sek - 40*40
  #pragma omp parallel
  {
    //Thread specific variables
    double start;
    double end;
    #pragma omp single
    cout << "Number of threads in use: " << omp_get_num_threads() << endl;

    #pragma omp for
    for (int i = 0; i < T_length; i++){
      start = omp_get_wtime();
      IsingModel is = IsingModel(L, T_array[i], 2); // n, temp, initmethod: (0)up, (1)down or (2)random
      //is.printMatrix();
      is.solve();
      end = omp_get_wtime();
      #pragma omp critical
      cout << "L: " << L << ".  Thread " << omp_get_thread_num() << " finished with: " <<  "T: " << T_array[i] << ". Time: " << end-start << "s" << endl;
      is.writeFile();

      //cout << "\n" << "----------------------------------------------" << endl;
    }
  }
  }



  /*
  IsingModel is = IsingModel(L, Ti, 2); // n, temp, initmethod: (0)up, (1)down or (2)random
  is.solve();
  is.printMatrix();
  TO DO:
  Run for L = 2, compare to analytical results.
  How does the number of accepted configurations behave as function of temperature T?
  Temp = 1 and 2.4 for L = 20 only

  Run simulation for multiple values of T
  Write Cv and chi to file for every T value.

  Spørsmål:
  Hvordan skrive verdi til 1 fil når det kjøret for flere tråder
  Skrive T, Cv, chi, M.

  */
  return 0;
}
