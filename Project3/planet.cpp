#include "planet.hpp"

using namespace std;
using namespace arma;

Planet::Planet(double m, double x, double y, double z, double vx, double vy, double vz, string thisName, int N_points){

  N = N_points;
  pos = vec(3*N); vel = vec(3*N);
  mass = m;
  name = thisName;
  pos[0] = x;
  pos[N] = y;
  pos[2*N] = z;
  vel[0] = vx;
  vel[N] = vy;
  vel[2*N] = vz;
  cout << "Initializing " << name << endl;

}


vec Planet::distanceOther(Planet& otherPlanet, int index, int N){
  //vec dr = pos - otherPlanet.pos;
  //return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));

  //compute distance to other planet for all three directions at a given time (=index)
  vec dis(3);
  dis[0] = pos[index] - otherPlanet.pos[index];
  dis[1] = pos[index+N] - otherPlanet.pos[index+N];
  dis[2] = pos[index+2*N] - otherPlanet.pos[index+2*N];

  //return  sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
  return dis;
}

vec Planet::gravitationalForce(Planet& otherPlanet, int index, int N){
  //double r = distanceOther(otherPlanet, index);
  //vec Fg = G_scale * mass * otherPlanet.mass / (r * r);
  //vec Fg =  otherPlanet.mass / (r * r); //use if TotalForceOnPlanet is used in solve

  vec r_vec =  - distanceOther(otherPlanet, index, N);
  double r = sqrt(dot(r_vec,r_vec));
  vec Fg = otherPlanet.mass * r_vec / (r*r*r); // !!!!!!! Rename Acceleration since not dividing by planet.mass?????????
  return Fg;
}
