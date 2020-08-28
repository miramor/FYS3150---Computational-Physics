
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

inline double exact(double x) { //u(x)
  return 1.0-(1-exp(-10))*x-exp(-10*x);
}

// d2_tilde = d2 - e1**2/d1
// g2_tilde = g2 - g1*e1/d1

// d2_tilde = d_i - e**2_i-1/
