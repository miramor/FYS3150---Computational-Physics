#include "SIRS.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>

using namespace std;

int main(int argc, char const *argv[]) {

  double S = 300;
  double I = 100;
  double a = 4;
  double b = 1;
  double c = .5;
  double t = 10;
  double dt = 1;

  //SIRS commonCold = SIRS(300., 100., 4., 1., 0.5, 1000 , 1);
  SIRS popA(S, I, a, b, c, t , dt);
  popA.solve("popA");
  return 0;
}
