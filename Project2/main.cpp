#include "jacobiAlgo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

double V0(double x);
double V1(double x);
double V2(double x);

int main(int argc, char const *argv[]) {

  JacobiEigenSolve jes;

  jes.Initialize(-1, 2, 3, V0); //parameter upperdiagonal, diagonal, matrix dimension, potential
  jes.Solve();
  //jes.PrintA();
  //jes.TestFindMaxEle();
  //jes.TestSolve();
  //jes.TestOrthogonality(); //writes out 1 if matrix is orthogonal else 0
  //jes.TestInitialize();

}

//Potential added along diagonal
double V0(double x){
  return 0;
}
double V1(double x){
  return x*x;
}

double V2(double x){
  double omega;
  omega = 0.25;
  return omega*omega*x*x + 1/x;
}
