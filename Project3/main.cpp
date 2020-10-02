#include "solver.hpp"
#include "planet.hpp"
#include <string.h>


int main (int argc, char const *argv[]){

  Planet sun(1.989e30, {0, 0, 0}, {0, 0, 0});
  Planet earth(5.972e24, {1, 0, 0}, {0, 1, 0});
  Planet jupyter(1.898e27, {2, 0, 0}, {0, 1, 0});
  return 0;
}
