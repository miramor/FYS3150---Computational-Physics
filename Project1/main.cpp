/*
Solving the one-dimensional Poisson equation with Dirchlet boundary conditions
Computing the solution with both a general algorithm and a specialised case, considering that the upper and lower diagonal are equal.
Algortihm includes forward and backward solution
*/

#include "algo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;


int main(int argc, char const *argv[]) {

  bool useSpecial = false;
  string filename = "CPUtime general";
  if(useSpecial){
    filename = "CPUtime special";
  }
  ofstream outfile;
  outfile.open(filename);

  for(int n=10; n<=10; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0); // a, b, c, n, x0, xn
    dSolv.Solve(useSpecial);
    dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Printtest();
    //cout << dSolv.solvetime << endl;
    //outfile << setprecision(6) << scientific;
    //outfile << n << ", " << dSolv.solvetime << endl;

  }
  outfile.close();


/*
  double a, b, c, x0, xn;
  a = -1.0; b = 2.0; c = -1.0; x0 = 0.0; xn = 1.0;
  double n = 100; //fjern
  DiffSolver dSolv;
  dSolv.Initialize(a, b, c, n, x0,  xn);
  dSolv.SolveLU(a, b, c);
  cout << "Time LUsolve " <<dSolv.solvetimeLU << endl;
*/


//Run SolveLU without Initialize?
//set precision
// loop for n?, read in max exponent using pow(10.0, i)
// string <<,<< instead string(comma)
//scientific notation when writing to csv file: WritetoFile()
// setiosflags(ios::showpoint | ios::uppercase) for scientific notation
//cpp plus for string formatation
//do we use LU decomposition?








  return 0;
}
