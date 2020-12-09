#include "SIRS.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <list>
#include <map>

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {

  double S = 300;
  double I = 100;
  double a = 4;
  //double b = 1;
  double c = .5;
  double e = 0.1;
  double d= 0.1;
  double dI = 0.1;
  double t_MC = 3;
  double t_RK4 = 20;
  double dt = 0.0025;
  int MC_cycles = 1000;

  list<double> b_val = {1, 2, 3, 4};
  map<double, double> b_totimeRK4;
  map<double, double> b_totimeMC;

  b_totimeRK4[1.] = 20.;
  b_totimeRK4[2.] = 20.;
  b_totimeRK4[3.] = 20.;
  b_totimeRK4[4.] = 20.;

  b_totimeMC[1.] = 200.;
  b_totimeMC[2.] = 200.;
  b_totimeMC[3.] = 200.;
  b_totimeMC[4.] = 200.;
  //list<double> b_val = {1};
  //list<double> b_val = {2};

  for(double b : b_val){
    //SIRS popMC(S, I, a, b, c, t_MC , MC_cycles);
    // SIRS pop(S, I, a, b, c, b_totimeMC[b], MC_cycles);
    // pop.solveMC("./Results/pop_" + to_string((int)b));

    SIRS popVD(S, I, a, b, c, b_totimeMC[b], MC_cycles, e, d, dI);
    popVD.solveMC("./Results/popVD_" + to_string((int)b));

    //SIRS popMC(S, I, a, b, c, t_RK4 , dt);
    //SIRS popRK4(S, I, a, b, c,  b_totimeRK4[b] , dt);
    //popRK4.solveRK4("./Results/pop_" + to_string((int)b));
  }


  return 0;
}
