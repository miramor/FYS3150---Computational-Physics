#include "jacobiAlgo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

int main(int argc, char const *argv[]) {

  JacobiEigenSolve jes;
  int n = 4;
  string filename = "results_num"
  jes.Initialize(-1, 2, n, "regular");
  jes.Solve();
  jes.PrintA();
  //jes.TestFindMaxEle();
  //jes.TestInitialize();

  jes.Write_Results(filename + to_string(4));

}
