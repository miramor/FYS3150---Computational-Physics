#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "jacobiAlgo.hpp"

using namespace std;
using namespace arma;


void JacobiEigenSolve::Initialize(double a_val, double b_val, int max_ite, int n_val){
  //Set class variables
  a = a_val; b = b_val; max_iterations = max_ite; n = n_val;
  max_iterations = (double) n * (double) n * (double) n;
  A = zeros<mat>(n,n);
  V = zeros<mat>(n,n);

  //Use this to in a test to compare
  /*vec X(n); X.fill(a_val);
  vec Y(n); Y.fill(b_val);
  A = toeplitz(X,Y);
  A = mat(n,n);
  */
  //Setup the A matrix.
  // Not finished, find a better way to fill it
  for(int d=0 ; d < n ; d++){
    A(d,d)=b_val;
    for(int j=d+1; j < n; j++){
      A(d,j) = a_val;
        A(j,d)=A(d,j);
    }
  }
  cout << A << endl;
}

// Returns a tuple of value of element and its index. (val, row, column) => f.eks (2.2, 1, 3)
tuple<double, int, int> JacobiEigenSolve::FindMaxEle(Mat<double> A){
  double max = 0.0;
  for (int i =0; i <n; i++){
    for (int j = i+1; j<n; j++){
      if ( fabs(A(i,j)) > max){
        max = fabs(A(i,j));
      }
    }
  }
  return {maks, i, j}

}

// finds the value of cos(theta) and sin(theta) and fills up V each time
void JacobiEigenSolve::Rotate(Mat<double> A, int row, int col){

}

 //Runs the rotation until we reached the max ite or reached the eps
void JacobiEigenSolve::Solve(){

}

// Prints A, used for checks
Mat<double> JacobiEigenSolve::PrintA(){
  cout << "Matrix A: " << A << endl;
}
