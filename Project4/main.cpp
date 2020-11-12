
#include <cstdio>
//#include <omp.h>
#include <cmath>
#include <iostream>
#include "isingModel.hpp"

using namespace std;


int main(int argc, char const *argv[]) {
  IsingModel is = IsingModel(20,2.4, 2); // n, temp, initmethod: (0)up, (1)down or (2)random
  //is.printMatrix();
  is.solve();






  /*
  TO DO:
  Run for L = 2, compare to analytical results.
  How does the number of accepted configurations behave as function of temperature T?
  Temp = 1 and 2.4 for L = 20 only

  Run simulation for multiple values of T
  Write Cv and chi to file for every T value.

  */
  return 0;
}
