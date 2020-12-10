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
  double t_MC = 2.;
  double t_RK4 = 2.;
  double dt = 0.0025;
  int MC_cycles = 2;

  double e = 0.015;
  double d = 0.012;
  double dI = 0.1;

  //list<double> b_val = {1, 2, 3, 4};
  map<double, double> b_totimeRK4;
  map<double, double> b_totimeMC;

  b_totimeRK4[1.] = 20.;
  b_totimeRK4[2.] = 20.;
  b_totimeRK4[3.] = 20.;
  b_totimeRK4[4.] = 20.;

  b_totimeMC[1.] = 20;
  b_totimeMC[2.] = 20.;
  b_totimeMC[3.] = 20.;
  b_totimeMC[4.] = 20.;
  list<double> b_val = {1,2,3,4};

  for(double b : b_val){
    //SIRS popMC(S, I, a, b, c, t_MC , MC_cycles);
    SIRS popMC(S, I, a, b, c, b_totimeMC[b]);
    //popMC.specMC(MC_cycles);
    //popMC.specMC_VD(MC_cycles, e, d, dI);
    //popMC.solveMC("./Results/pop_" + to_string((int)b));

    SIRS popRK4(S, I, a, b, c, b_totimeRK4[b]);
    //popRK4.specRK4(dt);
    popMC.specRK4_VD(dt, e, d, dI);
    popRK4.solveRK4("./Results/pop_" + to_string((int)b));
  }


  return 0;
}
