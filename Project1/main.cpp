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
  // First half of writes the time to compute the algorithm to text files
  // Second uses armadillo solve function (which uses LU?)
  // Can print out error and u(solution) before we write to file if we wanna check it for testing.
  bool useSpecial = true;
  int n_max = 1e7;
  int repeat = 8; // To avoid randomness, average out time by repeating calculation multiple times. Used 40 in report
  string filename = "CPUtime general";
  if(useSpecial){
    filename = "CPUtime special";
  }
  ofstream outfile;
  outfile.open(filename);


  for(int n=10; n<=n_max; n*= 10){
    DiffSolver dSolv;
    double totaltime = 0.0;

    for(int i = 0; i< repeat; i++){
      dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0, useSpecial); // Sincer we delete needed arrays in solve, need to initialize each time
      dSolv.Solve();
      totaltime += dSolv.solvetime;
    }

    dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Printtest();
    outfile << setprecision(7) << scientific;
    outfile << n << ", " << totaltime/(double)repeat <<endl;
    cout << "Finished Solve for n: " << n << endl;

  }
  outfile.close();*/

  //*****************************************************************************************************************

  // Run up to 1e4 for LU decomp, due to error in LU when to big matrices
  // Keep in mind if many repeats the process will take a long time for 1e4
  //ofstream outfile;

  outfile.open("CPUtime LU");

  for(int n=10; n<=1e4; n*= 10){
    DiffSolver dSolv;

    double totaltime = 0.0;

    for(int i = 0; i< repeat; i++){
      dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0, useSpecial); // a, b, c, n, x0, xn
      dSolv.SolveLU(-1.0, 2.0, -1.0);
      totaltime += dSolv.solvetimeLU;
    }
    outfile << setprecision(7) << scientific;
    outfile << n << ", " << totaltime/(double)repeat <<endl;
    cout << "Finished LU for n: " << n << endl;

  }
  outfile.close();

  return 0;
}
