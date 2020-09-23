#include "jacobiAlgo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

int main(int argc, char const *argv[]) {

  JacobiEigenSolve jes;
  int n = atoi(argv[1]);
  string filename = "results_num";
  string solution = argv[2];
  cout << solution << endl;
  jes.Initialize(-1, 2, n, solution);
  jes.Solve();
  jes.PrintA();
  //jes.TestFindMaxEle();
  //jes.TestInitialize();

  jes.Write_Results(filename + to_string(n), "numerical");

}
