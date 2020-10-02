#include "planet.hpp"

Planet::planet(vec pos){

}



Planet::distanceOther(){
  vec dr = this->pos - otherPlanet.pos;
  return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));
}

Planet::gravitationalForce(Planet otherPlanet){
  double r = distanceOther(Planet otherPlanet);
  vec Fg = this->G * this->mass * otherPlanet.mass / (r * r);
  return Fg
}
