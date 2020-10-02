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
    vec pos(3);
    vec vel(3);
    double kE;
    double pE;

  public:
    double planet(double mass, vec position, vec velo);
    vec distanceOther(Planet otherPlanet);
    vec gravitationForce(Planet otherPlanet);
    double kineticEnergy();
    double potentialEnergy();


}
