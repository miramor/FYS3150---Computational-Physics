
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
  int L = atoi(argv[1]);
  double Ti = stod(argv[2]);
  double Tf = stod(argv[3]);
  int n_T = atoi(argv[4]);
  long int num_cycles= stol(argv[5]);

  double dT = (Tf-Ti)/((double) n_T -1);
  double T_array [n_T];

  cout << "Simulating for temperatures: ";
  for (int i = 0; i < n_T; i++){
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

  //Asks user for num of threads to use
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
    //Each thread executes for a seperate temperature at a time
    for (int i = 0; i < n_T; i++){
      start = omp_get_wtime();
      double thread_seed = time(0) + omp_get_thread_num(); //unique seed for each object
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
  }
  return 0;
}
