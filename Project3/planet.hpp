#ifndef PLANET_HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>
using namespace arma;
using namespace std;

class Planet{
  private:
    vec acceleration;
    double G = 6.67408e-11; //Gravitational constant
    double pi = 2*acos(0.0);
    double G_scale = 4*pi*pi;
    int N;
    double c_sq = pow(63239.7263,2); //Speed of light squared in AU per year

  public:
    double l_merc; // ang moment for Mercury
    double mass;
    vec pos;
    vec vel;
    vec KEvec;
    vec PEvec;
    string name;

    Planet(double m, double x, double y, double z, double vx, double vy, double vz, string thisName, int N_points);
    vec distanceOther(Planet& otherPlanet, int index);
    vec gravitationalForce(Planet& otherPlanet, int index);

    // Methods for not storing data in planets
    Planet(double m, double x, double y, double z, double vx, double vy, double vz, string thisName, int N_points, int uselessNr);
    vec distanceOther_opt(Planet& otherPlanet, bool useCurr);
    vec gravitationalForce_opt(Planet& otherPlanet, bool useCurr);
    void update();

    Planet(double m, vec position, vec velocity);

    vec getKEvec();
    vec getPEvec();
};

#endif
