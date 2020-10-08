#ifndef PLANET_HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>
using arma::vec;

class Planet{
  private:
    vec acceleration;
    double kE;
    double pE;
    double G = 6.67408e-11; //Gravitational constant
    double pi = 2*acos(0.0);
    double G_scale = 4*pi*pi;
    int N = 3;


  public:

    double mass;
    vec pos;
    vec vel;
    std::string name;

    Planet(double m, vec position, vec velocity, std::string name);
    vec distanceOther(Planet& otherPlanet, int index, int N);
    vec gravitationalForce(Planet& otherPlanet, int index, int N);

    Planet(double m, vec position, vec velocity);
    double distanceOther(int N_val, const Planet& otherPlanet, int index);
    vec gravitationForce(int N_val, const Planet& otherPlanet, int index);

    double kineticEnergy();
    double potentialEnergy();
};

#endif
