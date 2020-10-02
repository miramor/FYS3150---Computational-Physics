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
    Planet(double m, vec position, vec velocity);
    double distanceOther(const Planet& otherPlanet);
    vec gravitationForce(const Planet& otherPlanet);
    double kineticEnergy();
    double potentialEnergy();

}

#endif
