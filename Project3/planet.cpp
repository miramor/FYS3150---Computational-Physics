#include "planet.hpp"

using namespace std;
using namespace arma;

Planet::Planet(double m, double x, double y, double z, double vx, double vy, double vz, string thisName, int N_points){
  N = N_points;
  KEvec = vec(N);
  PEvec = vec(N);
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


vec Planet::distanceOther(Planet& otherPlanet, int index){
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

vec Planet::gravitationalForce(Planet& otherPlanet, int index){
  //double r = distanceOther(otherPlanet, index);
  //vec Fg = G_scale * mass * otherPlanet.mass / (r * r);
  //vec Fg =  otherPlanet.mass / (r * r); //use if TotalForceOnPlanet is used in solve

  vec r_vec =  - distanceOther(otherPlanet, index);
  double Beta = 2;
  double r = norm(r_vec);
  double r_sqr = r*r;

  vec Fg;
  //if ((name == "Mercury" && otherPlanet.name = "Sun") || (name == "Sun" &&  otherPlanet.name = "Mercury")){
  if ( false == true){
    Fg = otherPlanet.mass * (r_vec/r) * (1 + 3*l_merc*l_merc/(r_sqr*c_sq)) / r_sqr;
    cout << "Rel force" << Fg << endl;
    cout << " force" << otherPlanet.mass * (r_vec/r) / (pow(r,Beta)) << endl;

    //cout << "This planet: "<< name << ", f: " << Fg << endl;
  }
  else{
    Fg = otherPlanet.mass * (r_vec/r) / (pow(r,Beta));
  }
  return Fg;
}


// double Planet::getKE(){
//   return KEvec;
// }
//
// double Planet::getPE(){
//   return PEvec;
// }
