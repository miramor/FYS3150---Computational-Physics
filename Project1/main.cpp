/*
Solving the one-dimensional Poisson equation with Dirchlet boundary conditions
Computing the solution with both a general algorithm and a specialised case, considering that the upper and lower diagonal are equal.
Algortihm includes forward and backward solution
*/

#include "algo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


int main(int argc, char const *argv[]) {

  ofstream outfile;
  string filename = "CPUtime general";
  //cout << "Filename: " <<filename << endl;
  //outfile.open(filename);

  for(int n=10; n<=10e3; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0);
    dSolv.Solve();
    dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Initialize(a_val = -1.0, b_val = 2.0, c_val = -1.0, n = 100, x_0 = 0.0, x_n = 1.0);
    //dSolv.Printtest();
    cout << dSolv.solvetime << endl;

    //string outline = "";
    //outline.append(to_string(n));
    //outline.append(", ");
    //outline.append(to_string(dSolv.solvetime));

    //outfile << outline << endl; // Add final line to end of text file
  }
  //outfile.close();


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

/* CPU time
 n=10: general: special:
 n=100: general: special:
 n=1000: general: special:


*/
