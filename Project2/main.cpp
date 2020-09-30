#include "jacobiAlgo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

double V0(double x, double omega);
double V1(double x, double omega);
double V2(double x, double omega);
void writeTimeResults(double omega);

int main(int argc, char const *argv[]) {


  JacobiEigenSolve jes;
  string filename;
  int n = atoi(argv[1]);
  string solution = argv[2];
  double rho_max = stod(argv[3]);
  double omega = stod(argv[4]);

  if (solution =="V2"){
    filename = "_n_" + to_string(n) + "_" + solution + "_w_" + to_string(omega) + ".txt";
  }
  else {
    filename = "_n_" + to_string(n) + "_" + solution + ".txt";
  }
  if (solution == "V0")
      jes.Initialize(-1, 2, n, 1, V0, omega); //parameter upperdiagonal, diagonal, matrix dimension, rhomax, potential
  else if (solution == "V1")
      jes.Initialize(-1, 2, n, rho_max, V1, omega);
  else if (solution == "V2")
      jes.Initialize(-1, 2, n, rho_max, V2, omega);

  //jes.Solve();

  //jes.Write_Results(filename, solution);
  //jes.PrintA();
  //jes.TestFindMaxEle();
  //jes.TestSolve();
  //jes.TestOrthogonality(); //writes out 1 if matrix is orthogonal else 0
  writeTimeResults(omega);
  //jes.TestInitialize();
  //writeTimeResults(omega);
}

//Potential added along diagonal
double V0(double x, double omega){
  return 0;
}
double V1(double x, double omega){
  return x*x;
}

double V2(double x, double omega){
  return omega*omega*x*x + 1/x;
}

void writeTimeResults(double omega){
  ofstream outfile;
  outfile.open("results/V0/TimeTable.csv");
  outfile << setprecision(4);

  for(int n = 10; n <= 200; n += 10){
    for (int j = 0; j < 10; j ++){
      JacobiEigenSolve jes;
      jes.Initialize(-1, 2, n, 1, V0, omega);
      jes.Solve();
      outfile << n << "," << jes.timeArma << "," << jes.timeClass << "," << jes.iterations << endl;
    }
  }
  return;
}
