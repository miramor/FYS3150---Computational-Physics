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


  public:

    double mass;
    vec pos;
    vec vel;
    vec KEvec;
    vec PEvec;
    string name;

    Planet(double m, double x, double y, double z, double vx, double vy, double vz, string thisName, int N_points);
    vec distanceOther(Planet& otherPlanet, int index);
    vec gravitationalForce(Planet& otherPlanet, int index);

    Planet(double m, vec position, vec velocity);

    vec getKEvec();
    vec getPEvec();
};

#endif
