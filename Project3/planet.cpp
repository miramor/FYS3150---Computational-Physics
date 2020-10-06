#include "planet.hpp"

using namespace std;
using namespace arma;

Planet::Planet(string name, double m, vector<double> position, vector<double> velocity){
  name = name;
  mass = m;
  pos = position;
  vel = velocity;
}



Planet::distanceOther(const Planet& otherPlanet){
  vector dr = pos - otherPlanet.pos;
  return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));
}

Planet::gravitationalForce(const Planet& otherPlanet){

  double r = distanceOther(otherPlanet);
  vector<double> Fg = G * mass * otherPlanet.mass / (r * r);
  return Fg
}
