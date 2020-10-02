#ifndef PLANET.HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>

using namespace std;
using arma::vec;

class Planet{
  private:
    double mass;
    vec pos(3);
    vec vel(3);
    double kE;
    double pE;

  public:
    planet(double mass, vec position, vec velo);
    distanceOther(Planet otherPlanet)
    GravitationForce()
    kineticEnergy()
    potentialEnergy()


}
