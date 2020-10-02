#ifndef PLANET.HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>

using std::double;
using arma::vec;

class Planet{
  private:
    double mass;
    vec<double> pos;
    vec<double> vel;
    vec<double> acceleration;
    double kE;
    double pE;
    double G = 6.67408e-11 //Gravitational constant

  public:
    double planet(double mass, vec<double> position, vec<double> velo);
    vec distanceOther(Planet otherPlanet);
    vec gravitationForce(Planet otherPlanet);
    double kineticEnergy();
    double potentialEnergy();

}

#endif
