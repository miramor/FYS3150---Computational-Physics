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


  public:
    Planet(string name, double m, vector<double> position, vector<double> velocity);
    double distanceOther(const Planet& otherPlanet);
    vector<double> gravitationForce(const Planet& otherPlanet);
    double kineticEnergy();
    double potentialEnergy();
};

#endif
