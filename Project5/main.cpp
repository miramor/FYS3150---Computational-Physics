#include "SIRS.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <list>

using namespace std;

int main(int argc, char const *argv[]) {

  double S = 300;
  double I = 100;
  double a = 4;
  //double b = 1;
  double c = .5;
  double t = 20;
  double dt = 0.001;

  list<double> b_val = {1, 2, 3, 4};
  //SIRS commonCold = SIRS(300., 100., 4., 1., 0.5, 1000 , 1);
  for(double b : b_val){
    SIRS pop(S, I, a, b, c, t , dt);
    pop.solve("./Results/pop_" + to_string((int)b));
  }


  return 0;
}
