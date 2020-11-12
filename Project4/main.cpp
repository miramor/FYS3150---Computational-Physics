
#include <cstdio>
//#include <omp.h>
#include <cmath>
#include <iostream>
#include "isingModel.hpp"

using namespace std;


int main(int argc, char const *argv[]) {
  double L = atoi(argv[1]);
  double Ti = stod(argv[2]);
  double Tf = stod(argv[3]);
  double dT = stod(argv[4]);
  IsingModel is = IsingModel(L, Ti, 2); // n, temp, initmethod: (0)up, (1)down or (2)random
  is.solve();


  /*
  for {double i = Ti; i <= Tf; i+=dT}{
    IsingModel is = IsingModel(L,i,2); // n, temp, initmethod: (0)up, (1)down or (2)random
    //is.printMatrix();
    is.solve();
  }
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
