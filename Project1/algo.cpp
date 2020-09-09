
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


// The relevant functions used in the project. g = f*hˆ2 and the exact solution
inline double f(double x) {
  return 100.0*exp(-10.0*x);
}

inline double exact(double x) { //u(x)
  return 1.0-(1-exp(-10))*x-exp(-10*x);
}

/*
 * Initialize the object with all the needed paramaters and stores them.
 *
 * Create new arrays with n places for a,b,c,g,x and error.
 * Make arrays with size n+1 for solution and u.
 * Compute g = f(x)*hˆ2 and fills arrays. Precaculates b[i] if we useSpecial
 * in hopes to make the solve algorithm faster.
 *
 * @params a_val : lower diagonal, b_val : middle diagonal, c_val : upper diagonal,
 * n_val : number of points, x_0_val : startpoint, x_n_val : endpoint,
 * useSpecial : choice of using specialised Thomas algorithm if a_val == c_val.
 */
void DiffSolver::Initialize(double a_val, double b_val, double c_val, int n_val, double x_0_val, double x_n_val, bool useSpecial){
  x_0 = x_0_val;
  x_n = x_n_val;
  n = n_val;
  m_useSpecial = useSpecial;
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
      if (m_useSpecial){
        b[i] = (i+1.0)/i;
      }
      else{
        b[i] = b_val;
      }
      a[i] = a_val;
      c[i] = c_val;
  }

  exact_[0] = exact(x_0);
  exact_[n] = exact(x_n);
}

/*
 * Solves using the Thomas algorithm, special or general.
 *
 * Using forward and backwards substitution it calculates the solution
 * given the type of method (special or general) which is taken from
 * the class variable m_useSpecial. Also measure the time and stores it
 * in the class as a variable.
 */
void DiffSolver::Solve(){
  clock_t start, finish, finish2; // declare start and final time
  start = clock();

  for(int i = 2; i < n; i++){
    if (m_useSpecial){
      //b[i] = (i+1.0)/i; XXXXXXXXXXXX
      g[i] = g[i] + g[i-1]/b[i-1];
    }
    else{
      b[i] = b[i] - a[i-1]*c[i-1]/b[i-1];
      g[i] = g[i] - a[i-1]*g[i-1]/b[i-1];
    }
  }

  u[0] = 0; //boundry condition
  u[n] = 0; //boundry condition
  u[n-1] = g[n-1]/b[n-1]; //boundry condition where u[n] = 0

  int i = n - 2;
  while (i>0){
    if (m_useSpecial){
      u[i] = (g[i] + u[i+1])/b[i]; //3FLOPSx(n-2)
    }
    else{
      u[i] = (g[i] - c[i]*u[i+1])/b[i]; //3FLOPSx(n-2)
    }
    i -= 1;
  }
  //Fill array with log of relative error:
  for(int i = 1; i<n; i++){
    double relError = fabs( (u[i] - exact_[i])/exact_[i] );
    error[i] = log10(relError);
  }

  finish = clock();
  solvetime = ( (finish - start)/(double)CLOCKS_PER_SEC );
  //cout <<"Time to solve: (s): "<<( (finish - start)/(double)CLOCKS_PER_SEC ) << endl;

  //Delete unneccesary arrays: a, b, c, g
  delete [] a; delete [] b; delete [] c; delete [] g;
}

/*
 * Simple function to print error (used for tests)
 */
void DiffSolver::PrintError(){
  for (int i = 1; i<n; i++){
    cout << setprecision(8) << error[i] << endl;
  }
}

/**
 * Writes the results from using the Solve function creates.
 *
 * Writes a comma seperated file containing relevant info.
 * Row format:  x_i, solution, exact_, log_rel_error
 * Then deletes the arrays containing x and error to free up memory no longer needed.
 */
void DiffSolver::WritetoFile(){
  // Write to file (CSV file) with 4 columns: [x_i, solution, exact_, log_rel_error]
  ofstream outfile;
  string filename ="ResultsComputation/ResultsG_nval=" ;
  if (m_useSpecial){
    filename = "ResultsComputation/ResultsS_nval=";
  }
  filename.append(to_string(n));
  outfile.open(filename);

  //Write the results into a comma seperated file to make use of pandas read_csv()
  for(int i = 1; i < n; i++){
    outfile << setprecision(7);
    outfile << x[i] << ", " << u[i] << ", " << exact_[i] << ", " << error[i] << endl;
  }
  outfile.close();

  //Delete no longer needed data, since its available in a file
  delete[] x; delete[] error;
}

/*
 * Simple function to print out the results and the exact solution to compare.
 * Used for verifying results when making small changes.
 */
void DiffSolver::Printtest(){
  // Simple test to verify result
  for(int i = 1; i < n; i+= 1){
    cout <<"Sol: " <<u[i] << "   Exact: " << exact_[i] <<endl;
  }
}

/*
 * Solve using LU - decomposition, using armadillos built in functions.
 *
 * Fill matrix with the 3 values given. Decompose the matrix in two tridiagonal
 * matrices L and U. Solves the equations L*g = y and U*y = solution.
 * Then writes these results to the ResultsComputation folder.
 * Also measure time taken and stores the value in class variable.
 *
 * @param  a = values in the lower diagonal, b = middle diagonal, c = upper diagonal
 */
void DiffSolver::SolveLU(double a_val, double b_val, double c_val){
  clock_t start, finish; // declare start and final time
  start = clock();

  int n_ = n-1; //Make sure dont change the class variable n
  mat A = zeros<mat>(n_,n_);
  // Set up arrays for the simple case
  vec g(n_);  vec x(n_); //Ax=g

  A(0,0) = b_val;  A(0,1) = c_val;  x(0) = h;  g(0) =  h_sq*f(x(0)); //Correction of cases where i+1 gives error in for loop
  x(n_-1) = x(0)+(n_-1)*h; g(n_-1) = h_sq*f(x(n_-1));
  for (int i = 1; i < n_-1; i++){
    x(i) = x(i-1)+h;
    g(i) = h_sq*f(x(i));
    A(i,i-1)  = a_val;
    A(i,i)    = b_val;
    A(i,i+1)  = c_val;
  }
  A(n_-1,n_-1) = b_val; A(n_-2,n_-1) = a_val; A(n_-1,n_-2) = c_val;
  //A.print("A ="); // Check if matrix is set up correctly
  mat L, U;
  lu(L,U,A); //find LU decomposition
  //(A-L*U).print("Test of LU decomposition");   //Check that A = LU, which means we should get 0 - matrix

  vec y = solve(L,g); // find y, Ly=g there y=Ux using forward substitution
  vec solution = solve(U,y); // find x, Ux=y using backward substitution
  finish = clock();

  //Writes to file to compare
  ofstream outfile;
  string filename ="ResultsComputation/ResultsLU_nval=";
  filename.append(to_string(n));
  outfile.open(filename);
  outfile << solution << endl;
  outfile.close();

  solvetimeLU = ( (finish - start)/(double)CLOCKS_PER_SEC );
}
