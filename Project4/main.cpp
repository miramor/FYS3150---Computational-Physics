
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "isingModel.hpp"

using namespace std;


int main(int argc, char const *argv[]) {
  int L = atoi(argv[1]);
  int Ti = (int)stod(argv[2])*100;
  int Tf = (int)stod(argv[3])*100;
  int dT = (int)(stod(argv[4])*100);

  ofstream Lfile;
  Lfile.open("Observables_" + to_string(L) + ".csv");
  Lfile <<  "T, <E>, <M>, Cv, chi" << endl;

  double start;
  double end;

  #pragma omp parallel
    {
      double start;
      double end;
      #pragma omp single
        cout << "Number of threads in use: " << omp_get_num_threads() << endl;

      #pragma omp for
      for (int i = Ti; i <= Tf; i+= dT){
        start = omp_get_wtime();
        IsingModel is = IsingModel(L, (double)i/100, 2); // n, temp, initmethod: (0)up, (1)down or (2)random
        //is.printMatrix();
        is.solve();
        end = omp_get_wtime();
        #pragma omp critical
          cout << "Thread " << omp_get_thread_num() << " finished with: " <<  "T: " << (double)i/100 << ". Time: " << end << "s" << endl;
          is.writeFile();

        //cout << "\n" << "----------------------------------------------" << endl;
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
