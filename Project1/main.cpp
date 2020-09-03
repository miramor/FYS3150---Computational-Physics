#include "algo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


int main(int argc, char const *argv[]) {
/*
  ofstream outfile;
  string filename = "CPUtime general";
  cout << "Filename: " <<filename << endl;
  outfile.open(filename);

  for(int n=10; n<=10e5; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0);
    dSolv.Solve();
    //dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Initialize(a_val = -1.0, b_val = 2.0, c_val = -1.0, n = 100, x_0 = 0.0, x_n = 1.0);
    //dSolv.Printtest();
    cout << dSolv.solvetime << endl;

    string outline = "";
    outline.append(to_string(n));
    outline.append(", ");
    outline.append(to_string(dSolv.solvetime));

    outfile << outline << endl; // Add final line to end of text file
  }
  outfile.close();

*/
  double a, b, c, x0, xn;
  a = -1.0; b = 2.0; c = -1.0; x0 = 0.0; xn = 1.0;
  double n = 10000; //fjern
  DiffSolver dSolv;
  dSolv.Initialize(a, b, c, n, x0,  xn);
  dSolv.Solve();
  dSolv.SolveLU(a, b, c);
  cout << "Time solve " <<dSolv.solvetime << endl;
  cout << "Time LUsolve " <<dSolv.solvetimeLU << endl;

  return 0;
}

/* CPU time
 n=10: general: special:
 n=100: general: special:
 n=1000: general: special:


*/
