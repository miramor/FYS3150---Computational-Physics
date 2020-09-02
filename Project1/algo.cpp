
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "algo.hpp"
#include "time.h" // you have to include the time.h header

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

void DiffSolver::Initialize(double a_val, double b_val, double c_val, int n_val, double x_0_val, double x_n_val){
  x_0 = x_0_val;
  x_n = x_n_val;
  n = n_val;
  a = new double[n];
  b = new double[n]; // Diagonal
  c = new double[n];
  g = new double[n]; // Right side

  x = new double[n];
  error = new double[n];

  u = new double[n+1];//Solution
  exact_ = new double[n+1];// Exact solution
  double h_val = (x_n-x_0)/(double)n;
  h = h_val;
  h_sq = pow( h_val , 2);

  for(int i = 1; i<n; i++){
      double x_i = x_0+i*h;
      x[i] = x_i;
      exact_[i] = exact(x_i);
      g[i] = f(x_i)*h_sq;
      b[i] = b_val;
      a[i] = a_val;
      c[i] = c_val;
  }

  exact_[0] = exact(x_0);
  exact_[n] = exact(x_n);
}

void DiffSolver::Solve(){

  clock_t start, finish, finish2; // declare start and final time
  start = clock();


  for(int i = 2; i < n; i++){
      if (a[0] != c[0]){
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


  //Fill array with log of relative error:
  for(int i = 1; i<n; i++){
    double relError = fabs( (u[i] - exact_[i])/exact_[i] );
    //error[i] = 1.4;
    error[i] = log10(relError);
  }

  }

  finish = clock();
  solvetime = ( (finish - start)/(double)CLOCKS_PER_SEC );
  cout <<"Time to solve: (s): "<<( (finish - start)/(double)CLOCKS_PER_SEC ) << endl;

  finish2 = clock();
  //Delete unneccesary vectors: a, b, c, g
  delete [] a; delete [] b; delete [] c; delete [] g;
  cout <<"Time with delete (s): "<< ( (finish2 - start)/(double)CLOCKS_PER_SEC ) << endl;
  //cout << *(a) << endl;
  //cout << *(a+1) << endl;
  //cout << *(b) << endl;
  //cout << *(b+1) << endl;
}

void DiffSolver::PrintError(){
  for (int i = 1; i<n; i++){
    cout << setprecision(8) << error[i] << endl;
  }
}

void DiffSolver::WritetoFile(){
  // Write to file (CSV file) with 4 columns: [x_i, solution, exact_, log_rel_error]
  ofstream outfile;
  string filename = "Results nval = ";
  filename.append(to_string(n));
  cout << "Filename: " <<filename << endl;
  outfile.open(filename);

  string comma = ", ";
  for(int i = 1; i < n; i++){
    string outline = "";
    outline.append(to_string(x[i]));
    outline.append(comma);
    outline.append(to_string(u[i]));
    outline.append(comma);
    outline.append(to_string(exact_[i]));
    outline.append(comma);
    outline.append(to_string(error[i]));
    outfile << outline << endl; // Add final line to end of text file
  }
  outfile.close();

  //Delete no longer needed data, since its in a csv file
  delete[] x; delete[] error;
}

void DiffSolver::Printtest(){
  for(int i = 1; i < n; i+= 1){
    cout <<"Sol: " <<u[i] << "   Exact: " << exact_[i] <<endl;
  }
}

void DiffSolver::SolveLU(double a_val, double b_val, double c_val){
  clock_t start, finish; // declare start and final time
  start = clock();
  n = n-1;
  mat A = zeros<mat>(n,n);
  // Set up arrays for the simple case
  vec g(n);  vec x(n); //Ax=g
  cout <<h<< endl;
  A(0,0) = b_val;  A(0,1) = c_val;  x(0) = h;  g(0) =  h_sq*f(x(0));
  x(n-1) = x(0)+(n-1)*h; g(n-1) = h_sq*f(x(n-1));
  for (int i = 1; i < n-1; i++){
    x(i) = x(i-1)+h;
    g(i) = h_sq*f(x(i));
    A(i,i-1)  = a_val;
    A(i,i)    = b_val;
    A(i,i+1)  = c_val;
  }
  A(n-1,n-1) = b_val; A(n-2,n-1) = a_val; A(n-1,n-2) = c_val;

    // solve Ax = g
  vec solution  = solve(A,g);
  //cout << solution << endl;
  finish = clock();
  solvetimeLU = ( (finish - start)/(double)CLOCKS_PER_SEC );
  


}
// d2_tilde = d2 - e1**2/d1
// g2_tilde = g2 - g1*e1/d1

// d2_tilde = d_i - e**2_i-1/
