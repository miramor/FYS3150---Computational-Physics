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
  int row, col;
  for (int i =0; i <n; i++){
    for (int j = i+1; j<n; j++){
      if ( fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        row = i;
        col = j;
      }
    }
  }
  return {max, row, col};

}

// finds the value of cos(theta) and sin(theta) and rotating matrix A fills up V each time
void JacobiEigenSolve::Rotate(Mat<double> A, int l, int k){
  /*
  Finding cos(theta) and sin(theta)
  */
  double cos_, sin_, tan_, tau;
  if (A(l,k) != 0.0){
    tau = (A(l,l)-A(k,k)) / (2.0*A(k,l));

    if (tau > 0){
      tan_ = 1.0 / (tau + sqrt(1.0 + tau*tau));
    }
    else{
      tan_ = -1.0 / (tau + sqrt(1.0 + tau*tau));
    }

    cos_ = 1.0 / sqrt(1.0 + tan_*tan_);
    sin_ = cos_ * tan_;
  }
  else{
    cos_ = 1.0;
    sin_ = 0.0;
  }

  /*
  //Rotating matrix A
  */
  double a_kk, a_ll, a_ik, a_il;
  a_kk = A(k,k);
  a_ll = A(l,l);

  //Computing new matrix elements with indices k and l
  A(k,k) = cos_*cos_*a_kk - 2.0*cos_*sin_*A(k,l) + sin_*sin_*a_ll;
  A(l,l) = sin_*sin_*a_kk + 2.0*cos_*sin_*A(k,l) + cos_*cos_*a_ll;
  A(k,l) = 0.0; //
  A(l,k) = 0.0;

  //Computing the new remaining elements
  for (int i = 0; i <n; i++){
    if (i != k && i != l){
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = cos_*a_ik - sin_*a_il;
      A(k,i) = A(i,k);
      A(i,l) = cos_*a_il - sin_*a_ik;
      A(l,i) = A(i,l);
    }
cout << n << endl;
    //Computing new eigenvectors
    double r_ik, r_il;
    r_ik = R(i,k);
    r_il = R(i,l);
    R(i,k) = cos_*r_ik - sin_*r_il;
    R(i,l) = cos_*r_il + sin_*r_ik;
  }
  return;
}



 //Runs the rotation until we reached the max ite or reached the eps
void JacobiEigenSolve::Solve(){
  double max_value;
  int row, col;
  int iterations = 0;
  auto[max_value, row, col] = FindMaxEle(A);

  while (max_value < eps && iterations < max_iterations ){
    Rotate(row, col);
    auto[max_value, row, col] = FindMaxEle(A);
    iterations ++;
  }

  return;
}

// Prints A, used for checks
Mat<double> JacobiEigenSolve::PrintA(){
  cout << "Matrix A: " << A << endl;
}
