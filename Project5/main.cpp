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

  double S = 300; //inital number of susceptibles
  double I = 100; //inital number of infected
  double a = 4; //rate of transmission
  double c = .5; //rate of immunity loss
  double dt = 0.0025; //time-step for RK4
  int MC_cycles = 1000; //number of Monte Carlo cycles

  double e = 0.009; //birth rate
  double d = 0.0075; //death reate
  double dI = 1; //death rate due to illness

//parameters for oscillating rate of transmission
  double A = 1.5;
  double A0 = 4;
  double PI = 4*atan(1);
  double frequency = 0.08*(2*PI);

  bool useSV;
  bool useVD;
  bool useV;
  bool useSTD;


  cout << "Would you like to enable seasonal variation for the rate of transmission (a)?\n 1. Yes \n 0. No " << endl;
  cin >> useSV;
  //cin >> useV >> useSeasVar >> useSV;


  cout << "Would you like to enable vital dynamics (deaths & births)?\n 1. Yes \n 0. No " << endl;
  cin >> useVD;

  if(useVD == false){
    cout << "Would you like to enable vaccines (susceptible -> recovered)?\n 1. Yes \n 0. No " << endl;
    cin >> useV;
  }

  if(useV == false && useVD == false){
    useSTD = true;
  }
  cout << "VitalDyn = " << useVD << ", vaccin = "  << useV << ", std = " << useSTD <<  endl;


// 0, 0, 0 -> using standard
// 1, 0, 1 -> using vaccines with seasonal variation
// 1, 1, * -> using vital dynamics, with seasonal
// 0, 0, 1 -> using vaccines without seasonal variation
// 0, 1, * -> vital dynamics, no seasonal
// 1, 0, 0 -> using standard, with seasonal

  map<double, double> b_totimeRK4;
  map<double, double> b_totimeMC;
  //t_end for each group
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

    SIRS popRK4(S, I, a, b, c, b_totimeRK4[b]);
    SIRS popMC(S, I, a, b, c, b_totimeMC[b]);

    if(useSTD){
      cout << "Spec standard" << endl;
      popRK4.specRK4(dt);
      popMC.specMC(MC_cycles);
    }

    if(useSV){
      cout << "Enable seasonal variation" << endl;
      popRK4.enableSeasVar();
      popMC.enableSeasVar();
    }

    if(useVD){
      cout << "Spec Vital Dyna" << endl;
      popRK4.specRK4_VD(dt, e, d, dI);
      popMC.specMC_VD(MC_cycles, e, d, dI);
    }

    else if(useV){
      cout << "Spec Vaccines" << endl;
      popRK4.specRK4(dt,100);
      popMC.specMC(MC_cycles, 100);
    }

    popMC.solveMC("./Results/pop_" + to_string((int)b));
    popRK4.solveRK4("./Results/pop_" + to_string((int)b));
  }


  return 0;
}
