#include "planet.hpp"

using namespace std;
using namespace arma;

Planet::Planet(double m, vec position, vec velocity, string thisName) : mass{m}, name{thisName}{

  pos = zeros<vec>(N*3);
  vel = zeros<vec>(N*3);
  /*
  for(int i = 0; i < 3; i++){
    pos[i] = position[i];
    vel[i] = velocity[i];
  }
  */

}



vec Planet::distanceOther(Planet& otherPlanet, int index, int N){
  //vec dr = pos - otherPlanet.pos;
  //return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));

  //compute distance to other planet for all three directions at a given time (=index)
  double r_x = pos[index] - otherPlanet.pos[index];
  double r_y = pos[index+N] - otherPlanet.pos[index+N];
  double r_z = pos[index+2*N] - otherPlanet.pos[index+2*N];
  vec dis{r_x, r_y, r_z};

  //return  sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
  return dis;
}

vec Planet::gravitationalForce(Planet& otherPlanet, int index, int N){
  //double r = distanceOther(otherPlanet, index);
  //vec Fg = G_scale * mass * otherPlanet.mass / (r * r);
  //vec Fg =  otherPlanet.mass / (r * r); //use if TotalForceOnPlanet is used in solve
  vec r =  distanceOther(otherPlanet, index, N);
  vec Fg = otherPlanet.mass / (r*r); // !!!!!!! Rename Acceleration since not dividing by planet.mass?????????
  return Fg;
}
