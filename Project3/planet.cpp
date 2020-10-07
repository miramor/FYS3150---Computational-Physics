#include "planet.hpp"

using namespace std;
using namespace arma;

Planet::Planet(string name, double m, vector<double> position, vector<double> velocity){
  name = name;
  mass = m;
  pos = position;
  vel = velocity;
}



double Planet::distanceOther(int N_val, const Planet& otherPlanet, int index){
  //vec dr = pos - otherPlanet.pos;
  //return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));

  //compute distance to other planet for all three directions at a given time (=index)
  double r_x = pos[index] - otherPlanet.pos[index];
  double r_y = pos[index+N_val] - otherPlanet.pos[index+N_val];
  double r_z = pos[index+2*N_val] - otherPlanet.pos[index+2*N_val];
  vec dis{r_x, r_y, r_z};
  //return  sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
  return dis;
}

vec Planet::gravitationalForce(int N_val, const Planet& otherPlanet, int index){
  //double r = distanceOther(otherPlanet, index);
  //vec Fg = G_scale * mass * otherPlanet.mass / (r * r);
  //vec Fg =  otherPlanet.mass / (r * r); //use if TotalForceOnPlanet is used in solve

  vec r =  distanceOther(N_val, otherPlanet, index);
  vec Fg = otherPlanet.mass / (r*r); // !!!!!!! Rename Acceleration since not dividing by planet.mass?????????
  return Fg
}
