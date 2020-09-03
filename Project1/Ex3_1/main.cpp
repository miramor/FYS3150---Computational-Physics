#include "derivative.hpp"
#include <iostream>
#include <cmath>

using namespace std;

double f(double x);



int main(int argc, char const *argv[]) {
  //Parametre
  double a = sqrt(2);
  vector<double> h_d = {1, 0.1, 0.01, 0.0001, 1e-6};
  vector<float> h_f = {1, 0.1, 0.01, 0.0001, 1e-6};

  derivativeAtPoint dvp;
  dvp.Initialize(h_d, h_f, a);
  dvp.FwEuler(f);
  dvp.BackAndFw(f);
  dvp.PrintResult();

  return 0;
}

double f(double x){
  return (1/tan(x));
}
