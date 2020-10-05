#include "planet.hpp"

using namespace std;
using namespace arma;

Planet::Planet(double m, vec position, vec velocity){
  mass = m;
  pos = position;
  vel = velocity;
}



Planet::distanceOther(const Planet& otherPlanet){
  vec dr = pos - otherPlanet.pos;
  return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));
}

Planet::gravitationalForce(const Planet& otherPlanet){

  double r = distanceOther(otherPlanet);
  //vec Fg = G_scale * mass * otherPlanet.mass / (r * r);
  vec Fg =  otherPlanet.mass / (r * r); //use if TotalForceOnPlanet is used in solve
  return Fg
}
