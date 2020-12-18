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
  int MC_cycles = 1000;

  double e = 0.009;
  double d = 0.0075;
  double dI = 5;

  double A = 1.5; double A0 = 4; double frequency = 0.08*(2*PI);

  map<double, double> b_totimeRK4;
  map<double, double> b_totimeMC;

  b_totimeRK4[1.] = 20;
  b_totimeRK4[2.] = 20;
  b_totimeRK4[3.] = 50;
  b_totimeRK4[4.] = 30;

  b_totimeMC[1.] = 20;
  b_totimeMC[2.] = 20;
  b_totimeMC[3.] = 50;
  b_totimeMC[4.] = 30;

  list<double> b_val = {1,2,3,4};

  for(double b : b_val){
    //SIRS popMC(S, I, a, b, c, t_MC , MC_cycles);
    SIRS popMC(S, I, a, b, c, b_totimeMC[b]);
    //popMC.specMC(MC_cycles, 100);
    popMC.specMC(MC_cycles);
    //popMC.enableSeasVar();
    //popMC.specMC_VD(MC_cycles, e, d, dI);
    popMC.solveMC("./Results/pop_" + to_string((int)b));

    SIRS popRK4(S, I, a, b, c, b_totimeRK4[b]);
    //popRK4.specRK4(dt,100);
    popRK4.specRK4(dt);
    //popRK4.enableSeasVar();
    //popRK4.specRK4_VD(dt, e, d, dI);
    popRK4.solveRK4("./Results/pop_" + to_string((int)b));
  }


  return 0;
}
