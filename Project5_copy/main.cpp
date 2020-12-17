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
  int MC_cycles = 500;

  double e = 0.009;
  double d = 0.0075;

  map<double, double> b_totimeRK4;
  map<double, double> b_totimeMC;
  double time = 30.;
  b_totimeRK4[1.] = time;
  b_totimeRK4[2.] = time;
  b_totimeRK4[3.] = time;
  b_totimeRK4[4.] = time;

  b_totimeMC[1.] = time;
  b_totimeMC[2.] = time;
  b_totimeMC[3.] = time;
  b_totimeMC[4.] = time;

  list<double> b_val = {1,2,3,4};

  ofstream ofile;
  ofile.open("dI_data.csv");

  double dI = 1;
  double deaths;
  double maxF = 600;
  for (double f = 0; f <= maxF; f += 5){
    cout << "*************************************" << endl;
    cout << "f = " << f << endl;
    ofile << f << ", ";
    for(double b : b_val){
      //SIRS popMC(S, I, a, b, c, t_MC , MC_cycles);
      SIRS popMC(S, I, a, b, c, b_totimeMC[b]);
      //popMC.specMC(MC_cycles, f);
      //popMC.specMC(MC_cycles);
      //popMC.specMC_VD(MC_cycles, e, d, dI);
      popMC.specMC_VD_VAC(MC_cycles, e, d, dI, f);
      popMC.solveMC("./Results/pop_" + to_string((int)b));
      //deaths = popMC.deadDis_final;
      //cout << b << " val -> deaths: " << popMC.deadDis_final << endl;
      //cout <<"S, I , R: " <<  popMC.S_final << ", " << popMC.I_final << ", " << popMC.R_final << endl;
      ofile << popMC.deadDis_final << ", ";
      //ofile << popMC.S_final << ", " << popMC.I_final << ", " << popMC.R_final << ", ";

      //SIRS popRK4(S, I, a, b, c, b_totimeRK4[b]);
      //popRK4.specRK4(dt,75);
      //popRK4.specRK4(dt);
      //popRK4.specRK4_VD(dt, e, d, dI);
      //popRK4.solveRK4("./Results/pop_" + to_string((int)b));
    }
    ofile << endl;
  }



  return 0;
}
