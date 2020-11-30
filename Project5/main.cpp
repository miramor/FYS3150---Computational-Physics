#include "SIRS.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>

using namespace std;

int main(int argc, char const *argv[]) {

  SIRS commonCold = SIRS(300., 100., 4., 1., 0.5, 1000 , 1);
  commonCold.solve("popA.csv");
  return 0;
}
