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
  a = a_val; b = b_val; max_iterations = max_ite; n = n_val

  //Use this to in a test to compare
  /*vec X(n); X.fill(a_val);
  vec Y(n); Y.fill(b_val);
  A = toeplitz(X,Y);
  A = mat(n,n);
  */
  //Setup the A matrix.
  // Not finished, find a better way to fill it
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if( i == j){
        A[i][j] = b_val ;
      }
      if( i == j+1);
      }
    }
  }

}

// Returns a tuple of value of element and its index. (val, row, column) => f.eks (2.2, 1, 3)
tuple<double, int, int> FindMaxEle(Mat<double> A){

}

// finds the value of cos(theta) and sin(theta) and fills up V each time
void Rotate(Mat<double> A, int row, int col){

}

 //Runs the rotation until we reached the max ite or reached the eps
void Solve(){

}

// Prints A, used for checks
Mat<double> PrintA(){
  cout << "Matrix A: " << A << endl;
}
