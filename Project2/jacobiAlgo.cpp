#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "jacobiAlgo.hpp"
#include <cassert>

using namespace std;
using namespace arma;

/*
inline vec<double> f(int i, double h) {
}
*/

void JacobiEigenSolve::Initialize(double a_val, double b_val, int n_val){
  //Set class variables
  n = n_val; //Points
  h = 1.0/(n+1);  // Step size (x-x0)/N = (x-x0)/(n+1)
  a = a_val/(h*h); b = b_val/(h*h);
  max_iterations = (double) n * (double) n * (double) n;
  //max_iterations = 100;
  A = zeros<mat>(n,n);
  R.eye(n,n);


  for (int i=0; i<n; i++){
    if (i<n-1){
      A(i,i) = b;
      A(i,i+1) = a;
      A(i+1,i) = a;
    } else{
      A(i,i) = b;
    }
  }

  A_test = repmat(A, 1, 1); //Make a copy of A to be used in tests
  vec eig = eig_sym(A);
  //cout << "Fasit eigenvalues: \n" << eig_sym(A) << endl;
  //cout << "Fasit eigenvectors: \n " << eigs_gen(A, n) << endl;
  //cout << A << endl;

  // Opg d - potensial. Add potential on diagonal elements.
  bool usePotential = false; //implement this later
  if(usePotential){
    int N = n+1;
    int p0 = 0;
    int pmax = 10;
    h = (pmax - p0)/N;
    for(int i = 1; i < N; i++){
      A(i-1,i-1) += (p0 + i*h); // remember to change to squared
    }
  }
  return;
}

// Returns a tuple of value of element and its index. (val, row, column) => f.eks (2.2, 1, 3)
void JacobiEigenSolve::FindMaxEle(double& max, int& row, int& col){
  max = 0.0;
  for (int i =0; i <n; i++){
    for (int j = i+1; j<n; j++){
      if ( fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        row = i;
        col = j;
      }
    }
  }
  return;
  //return {max, row, col};

}

// finds the value of cos(theta) and sin(theta) and rotating matrix A fills up V each time
void JacobiEigenSolve::Rotate(int l, int k){
  /*
  Finding cos(theta) and sin(theta)
  */
  double cos_, sin_, tan_, tau;
  //cout << "I rotate: \n"<< A <<"\n"<< A << endl;
  //cout << "l: " <<l <<  " k: " << k << endl;;
  if (A(l,k) != 0.0){
    //cout << "k " << k<< endl;
    tau = (A(l,l)-A(k,k)) / (2.0*A(k,l));
    //cout << "tau " << tau<< endl;

    if (tau > 0){
      tan_ = 1.0 / (tau + sqrt(1.0 + tau*tau));
    }
    else{
      tan_ = -1.0 / (-tau + sqrt(1.0 + tau*tau));
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
      A(i,l) = cos_*a_il + sin_*a_ik;
      A(l,i) = A(i,l);
    }
//cout << n << endl;
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
  int iterations = 0; int row; int col;
  double max_val;
  FindMaxEle(max_val, row, col);
  //cout << A << endl;

  while (max_val > eps || iterations < max_iterations ){
    Rotate(row, col);
    FindMaxEle(max_val, row, col);
    //cout <<  "Kolonne" <<col_ << "row: "<< row_ << endl;
    iterations ++;
  }
  A.clean(eps); //Remove elements smaller than eps.
  return;
}

// Prints A, used for checks
void JacobiEigenSolve::PrintA(){
  cout << "Matrix A: \n" << A << endl;
  cout << "Matrix R: \n" << R << endl;
  return;
}


//Tests 2-3 times if method finds correct value and postion
void JacobiEigenSolve::TestFindMaxEle(){
  int row; int col;
  double max_val;
  arma_rng::set_seed_random();
  A = mat(4 ,4, arma::fill::randn);
  A= trimatu(A);//make matrix upper triangular
  FindMaxEle(max_val, row, col);
  //cout << "Find max value of this matrix:\n" << A << endl;

  //To make use of index_max remove all diag and take abs value of all elements
  A = trimatu(A);//make matrix upper triangular since only search upper
  A.diag().zeros();
  A = abs(A);

  cout << A << endl;
  uvec s = ind2sub( size(A), A.index_max() );

  //cout << "Armadillo:  " << "Row:" << s(0) << " Col: " << s(1) << endl;
  //cout << "Classfunc:  " << "Row:" << row << " Col: " << col << endl;
  if( row == s(0) && col == s(1)){
    cout << "Test Succesful! Correct index for max element" << endl;
  }
  else{
    cout << "Col and Row does not match: (Arma vs Classfunc)\n" << "Col: " << s(0) << ", " << col << "\nRow: " << s(1) << ", " << row << endl;
  }

  return;
}

//Tests if matrix is set up correctly initially
void JacobiEigenSolve::TestInitialize(){
  //Use this to in a test to compare
  vec X(n); X.fill(0);
  X(0) = 2; X(1) = -1;
  mat B = toeplitz(X);

  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if( fabs(A_test(i,j) - B(i,j)) > 1e-8 ){
        cout << "Error: A(i,j) - B(i,j) = " << A_test(i,j) - B(i,j) << endl;
      }
    }
  }
  return;
}

void JacobiEigenSolve::TestSolve(){
  // Sort egenverdiene
  // Sjekk om egenverdiene er riktig etter vi har kjørt Solve()
  // Enten kan vi oppretee nytt objekt her eller anta at vi har kjørt løsning
  // Ved å opprette nytt betyr det at vi kan teste når som helst.
  JacobiEigenSolve jes;
  vec sortedEign =  sort(A.diag());
  //cout << sortedEign << endl;
  //cout << eig << endl;
  for(int i = 0; i < eig.size(); i++){
    if( fabs(eig(i) - sortedEign(i) ) > eps ){
      cout << "Error the eig values does not match index: " << i << "  Exact" << eig(i) << " Solved: "<< sortedEign(i) << endl;
    }
  }
  cout << "Succesful test. Correct eigenvalues" << endl;
}
