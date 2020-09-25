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

  string filename = "_n_" + string(argv[1]) + "_" + string(argv[2]) + ".txt";
  string solution = argv[2];
  int n = atoi(argv[1]);
  int rho_max = 10;

  if (solution == "V0")
      jes.Initialize(-1, 2, n, rho_max, V0); //parameter upperdiagonal, diagonal, matrix dimension, rhomax, potential
  else if (solution == "V1")
      jes.Initialize(-1, 2, n, rho_max, V1);
  else if (solution == "V2")
      jes.Initialize(-1, 2, n, rho_max, V2);

  jes.Solve();
  jes.Write_Results(filename, solution);
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
