#ifndef PLANET_HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>

using arma::vec;
using namespace std;

class Planet{
  private:
    string name;
    double mass;
    vector<double> pos;
    vector<double> vel;
    vector<double> acceleration;
    double kE;
    double pE;
    double G = 6.67408e-11; //Gravitational constant
    double pi = 2*acos(0.0);
    double G_scale = 4*pi*pi;


  public:
    Planet(double m, vec position, vec velocity);
    double distanceOther(int N_val, const Planet& otherPlanet, int index);
    vec gravitationForce(int N_val, const Planet& otherPlanet, int index);
    double kineticEnergy();
    double potentialEnergy();
};

#endif
