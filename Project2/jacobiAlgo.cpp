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

void JacobiEigenSolve::Write_Results(string filename, string solution){
  ofstream ofile;
  ofile.open(filename);

  if (solution == "analytical"){
    //Define analytical eval and evec
    int eval = 1;
    int evec = 1;
  }
  else{
    eval = max(A);
    evec = R;
  }
  for (int i = 0; i < n; i ++)
      ofile << x(i) << " ";
  ofile << endl << eval;

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      ofile << evec(i,j) << " ";
    }
    ofile << endl;
  }
}


void JacobiEigenSolve::Initialize(double a_val, double b_val, int n_val, string init){
  //Set class variables
  n = n_val; //Points
  h = 1.0/(n+1);  // Step size (x-x0)/N = (x-x0)/(n+1)
  a = a_val/(h*h); b = b_val/(h*h);
  max_iterations = (double) n * (double) n * (double) n;
  //max_iterations = 100;
  x = zeros<rowvec>(n);
  A = zeros<mat>(n,n);
  R.eye(n,n);


  for (int i = 0; i < n-1; i++){
      x(i) = i*h;
      A(i,i) = b;
      A(i,i+1) = a;
      A(i+1,i) = a;
  }
  A(n-1,n-1) = b;
  x(n-1) = (n-1) * h;


  A_test = repmat(A, 1, 1);

  // Opg c - potensial. Add potential on diagonal elements.
  if(init == "potential"){
    int p0 = 0;
    int pmax = 10;
    h = (pmax - p0)/(n+1);
    for(int i = 0; i < n; i++){
      A(i,i) = A(i,i) + x[i]*x[i];
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
  arma_rng::set_seed_random();
  A = mat(4,4, arma::fill::randu);
  cout << "Find max value of this matrix:\n" << A << endl;
  //A = rand
  //Mat<double> C =
  //int i = index_max(A_test);
  //int j = index_min(A_test);
/*
  while(i == j){
    i = 2;
  }
  //index_max
  return;*/
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
  // Sjekk om egenverdiene er riktig.
  for (int i =0; i <n; i++){
    for (int j = i+1; j<n; j++){
      if( i != j){
        assert(abs(A(i,j)) < eps );
      }
    }
  }
  // Cross product  a x b  = null_vektor hvis de er paralellel.
  // vector.clean(eps)
  // Sjekk at alle elementene er 0. feks == vec zeros(n)
}
