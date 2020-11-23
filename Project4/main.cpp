
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "isingModel.hpp"
#include <time.h>

using namespace std;


int main(int argc, char const *argv[]) {
  //double Ti = stod(argv[2]);
  //double Tf = stod(argv[3]);
  //double dT = stod(argv[4]);

  //int L = atoi(argv[1]);
  IsingModel is = IsingModel(20, 1, 2); // n, temp, initmethod: (0)up, (1)down or (2)random
  // cv/beta = varians histo
  //is.printMatrix();
  is.solve();
  //cout << setprecision(6) << "Sigma: " << is.getSigma() << endl;

 /*
  int numThreadx8 = 1;
  double Ti = 2;
  double Tf = 2.35;
  double dT = (Tf-Ti)/(8*numThreadx8-1);
  cout << dT << endl;
  int T_length = 8*numThreadx8;
  double T_array [T_length];

  for (int i = 0; i < T_length; i++){
    T_array[i] = Ti + dT*i;
    cout << T_array[i] << ", ";
  }
  cout << " with a dT = " << dT << endl;

  ofstream Lfile;
  Lfile.open("Observables_" + to_string(L) + ".csv");
  Lfile <<  "T, <E>, <M>, Cv, chi" << endl;

  int threads_available;
  #pragma omp parallel
  {
    #pragma omp single
    {
      threads_available = omp_get_num_threads();
      cout << "You have " << threads_available << " threads available, how many do you want to use? (max " << threads_available << ")" << endl;
    }
  }

  int num_threads;
  cin >> num_threads;

  if(num_threads > threads_available || num_threads< 1){
    cout << num_threads << " is not possible. Automatically changed to max_threads = " << threads_available << endl;
    num_threads = threads_available;
  }

  omp_set_num_threads(num_threads);
  #pragma omp parallel
  {
    //Thread specific variables
    double start;
    double end;
    #pragma omp single
    {
      cout << "Number of threads in use: " << omp_get_num_threads() << endl;
      cout << "----------------------------------------------" << endl;
    }
    #pragma omp for
    for (int i = 0; i < n_T; i++){
      start = omp_get_wtime();
      double thread_seed = time(0) + omp_get_thread_num();
      IsingModel is = IsingModel(L, T_array[i], 2, num_cycles, thread_seed); // n, temp, initmethod: (0)up, (1)down or (2)random
      //is.printMatrix();
      is.solve();
      end = omp_get_wtime();
      #pragma omp critical
      {
        cout << "Thread " << omp_get_thread_num() << " finished with " <<  "T=" << T_array[i] << ". Time: " << end-start << "s" << endl;
        is.printValues();
        cout << "----------------------------------------------" << endl;
        is.writeFile();
      }
    }
  }*/




  /*
  IsingModel is = IsingModel(L, Ti, 2); // n, temp, initmethod: (0)up, (1)down or (2)random
  is.solve();
  is.printMatrix();
  TO DO:
  Run for L = 2, compare to analytical results.
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
