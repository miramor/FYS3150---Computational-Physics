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
  double t_MC = 20;
  double t_RK4 = 20;
  double dt = 0.0025;

  string method = "MC";
  list<double> b_val = {1, 2, 3, 4};
  map<double, double> b_totimeRK4;
  map<double, double> b_totimeMC;

  b_totimeRK4[1.] = 20.;
  b_totimeRK4[2.] = 20.;
  b_totimeRK4[3.] = 20.;
  b_totimeRK4[4.] = 20.;

  b_totimeMC[1.] = 300.;
  b_totimeMC[2.] = 300.;
  b_totimeMC[3.] = 300.;
  b_totimeMC[4.] = 15.;
  //list<double> b_val = {1};

  for(double b : b_val){
    //SIRS popMC(S, I, a, b, c, t_MC , dt);
    SIRS popMC(S, I, a, b, c, b_totimeMC[b] , dt);
    popMC.solve("./Results/pop_" + to_string((int)b), "MC");

    //SIRS popMC(S, I, a, b, c, t_RK4 , dt);
    SIRS popRK4(S, I, a, b, c,  b_totimeRK4[b] , dt);
    popRK4.solve("./Results/pop_" + to_string((int)b), "RK4");
  }


  return 0;
}
