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
  // First choose to use special or general algorithm
  // First half of writes the time to compute the algorith to text files
  // Second uses armadillo solve function (which uses LU?)
  // Can print out error and u(solution) before we write to file if we wanna check it for testing.
  bool useSpecial = false;
  int n_max = 1e6;
  int repeat = 8; // To avoid randomness, average out time by repeating calculation multiple times
  string filename = "CPUtime general";
  if(useSpecial){
    filename = "CPUtime special";
  }
  ofstream outfile;
  outfile.open(filename);

  for(int n=10; n<=n_max; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0, useSpecial); // a, b, c, n, x0, xn

    double totaltime = 0.0;
    for(int i = 0; i< repeat; i++){
      dSolv.Solve();
      totaltime += dSolv.solvetime;
    }


    dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Printtest();
    outfile << setprecision(6) << scientific;
    outfile << n << ", " << totaltime/(double)repeat <<endl;
    cout << "Finished for n: " << n << endl;

  }
  outfile.close();


  // Run up to 1e4 for LU decomp, due to error in LU when to big matrices
  ofstream outfile2;
  outfile2.open("CPUtime LU");

  for(int n=10; n<=1e4; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0, useSpecial); // a, b, c, n, x0, xn
    dSolv.SolveLU(-1.0, 2.0, -1.0);
    dSolv.WritetoFile();
    outfile2 << setprecision(6) << scientific;
    outfile2 << n << ", " << dSolv.solvetime <<endl;
    cout << "Finished LU for n: " << n << endl;

  }
  outfile.close();

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
