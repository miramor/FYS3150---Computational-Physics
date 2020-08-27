
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;


// vector v = [v1,v2,v3..., v_n]
// b_itilde = f(x)*h**2 = -v_i+1 + 2*v_i - v_i-1

inline double f(double x) {
  return 100.0*exp(-10.0*x);
}

inline double exact(double x) {
  return 1.0-(1-exp(-10))*x-exp(-10*x);
}
