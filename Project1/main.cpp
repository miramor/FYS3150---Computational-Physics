<<<<<<< HEAD
=======
<<<<<<< HEAD
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
=======
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb
/*
Solving the one-dimensional Poisson equation with Dirchlet boundary conditions
Computing the solution with both a general algorithm and a specialised case, considering that the upper and lower diagonal are equal.
Algortihm includes forward and backward solution
*/

#include "algo.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
<<<<<<< HEAD
#include <iomanip>
=======
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb

using namespace std;


int main(int argc, char const *argv[]) {

<<<<<<< HEAD
  bool useSpecial = false;
  string filename = "CPUtime general";
  if(useSpecial){
    filename = "CPUtime special";
  }
  ofstream outfile;
  outfile.open(filename);

  for(int n=10; n<=10; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0);
    dSolv.Solve(useSpecial);
    dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Initialize(a_val = -1.0, b_val = 2.0, c_val = -1.0, n = 100, x_0 = 0.0, x_n = 1.0);
    //dSolv.Printtest();
    //cout << dSolv.solvetime << endl;
    outfile << setprecision(6) << scientific;
    outfile << n << ", " << dSolv.solvetime << endl;
=======
  ofstream outfile;
  string filename = "CPUtime general";
  //cout << "Filename: " <<filename << endl;
  //outfile.open(filename);

  for(int n=10; n<=10e5; n*= 10){
    DiffSolver dSolv;
    dSolv.Initialize(-1.0, 2.0, -1.0, n,  0.0, 1.0);
    dSolv.Solve();
    //dSolv.WritetoFile();
    //dSolv.PrintError();
    //dSolv.Initialize(a_val = -1.0, b_val = 2.0, c_val = -1.0, n = 100, x_0 = 0.0, x_n = 1.0);
    //dSolv.Printtest();
    cout << dSolv.solvetime << endl;

>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb
    //string outline = "";
    //outline.append(to_string(n));
    //outline.append(", ");
    //outline.append(to_string(dSolv.solvetime));
<<<<<<< HEAD
    //outfile << outline << endl; // Add final line to end of text file
  }
  outfile.close();
=======

    //outfile << outline << endl; // Add final line to end of text file
  }
  //outfile.close();
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb


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
<<<<<<< HEAD
=======
>>>>>>> 54bab2448c3dffbb49242067a465c94e0886e423
>>>>>>> 16a87e405674f650190568fbb3a0e84d915a8fdb
