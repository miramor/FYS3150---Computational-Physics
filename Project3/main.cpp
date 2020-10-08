#include "solver.hpp"
#include "planet.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]) {

  Planet sun(1.989e30, {0, 0, 0}, {0, 0, 0}, "sun");
  Planet earth(5.972e24, {1, 0, 0}, {0, 1, 0}, "earth");

  vector<Planet> planets;
  planets.push_back(sun);
  planets.push_back(earth);

  Solver solv(planets,5, 10);
  solv.VelocityVerlet();


  //Planet jupyter =planet(1.898E27 kg, vec p(3) = {2, 0, 0}, vec v(3) = {0, 1, 0});

  return 0;
}
