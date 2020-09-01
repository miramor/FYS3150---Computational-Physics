
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "algo.hpp"

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

void DiffSolver::Initialize(double a, double b, double c, int n, double x_0, double x_n){
  a = new double[n];
  b = new double[n]; // Diagonal
  c = new double[n];
  g = new double[n]; // Right side

  u = new double[n+1];//Solution
  exact_ = new double[n+1];// Exact solution

  h_sq = pow( (x_n-x_0)/n , 2);

  for(int i = 1; i<n; i++){
      double x_i = x_0+i*h;
      exact_[i] = exact(x_i);
      g[i] = f(x_i)*h_sq;
      b[i] = 2.0;
      a[i] = -1.0;
      c[i] = -1.0;
  }

  exact_[0] = exact(x_0);
  exact_[n] = exact(x_n);
}

void DiffSolver::Solve(){
  for(int i = 2; i < n; i++){
      if (a[0] == c[0]){
        b[i] = b[i] - 1/b[i-1];
      }
      else{
        b[i] = b[i] - a[i-1]*c[i-1]/b[i-1];
      }
      g[i] = g[i] - a[i-1]*g[i-1]/b[i-1];
  }

  u[0] = 0; //boundry condition
  u[n] = 0; //boundry condition
  u[n-1] = g[n-1]/b[n-1]; //boundry condition where u[n] = 0

  int i = n - 2;
  while (i>0){
      u[i] = (g[i] - c[i]*u[i+1])/b[i];
      i -= 1;
  }
}

void DiffSolver::PrintError(){
  for (int i = 0, i<n, i++){
    cout << abs(u[i]-exact_[i]/exact_[i]) << endl;
  }
}

void DiffSolver::WritetoFile(string filename){
}


// d2_tilde = d2 - e1**2/d1
// g2_tilde = g2 - g1*e1/d1

// d2_tilde = d_i - e**2_i-1/
