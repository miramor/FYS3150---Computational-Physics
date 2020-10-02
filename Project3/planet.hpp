#ifndef PLANET_HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>

using arma::vec;

class Planet{
  private:
    double mass;
    vec pos;
    vec vel;
    vec acceleration;
    double kE;
    double pE;
    double G = 6.67408e-11; //Gravitational constant

  public:
    double planet(double mass, vec position, vec velo);
    vec distanceOther(Planet otherPlanet);
    vec gravitationForce(Planet otherPlanet);
    double kineticEnergy();
    double potentialEnergy();

}

#endif
