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

//Different Initializing to not store all values
Planet::Planet(double m, double x, double y, double z, double vx, double vy, double vz, string thisName, int N_points, int uselessNr){
  N = N_points;
  int timpeptsNeeded = 2;
  pos = vec(timpeptsNeeded*3, fill::zeros); vel = vec(timpeptsNeeded*3, fill::zeros);
  mass = m;
  name = thisName;
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  vel[0] = vx;
  vel[1] = vy;
  vel[2] = vz;
  cout << "Init_opt " << name << endl;
}

vec Planet::distanceOther_opt(Planet& otherPlanet, bool useCurr){ // if not curr, then distance previous
  int i;
  if(useCurr == true){
    i = 3;
  }
  else{
    i = 0;
  }
  vec dis(3);
  dis[0] = pos[i] - otherPlanet.pos[i];
  dis[1] = pos[i+1] - otherPlanet.pos[i+1];
  dis[2] = pos[i+2] - otherPlanet.pos[i+2];

  return dis;
}

vec Planet::gravitationalForce_opt(Planet& otherPlanet, bool useCurr){
  vec r_vec =  - distanceOther_opt(otherPlanet, useCurr);
  double r = norm(r_vec);
  vec Fg;

  if (true == true){
    Fg = otherPlanet.mass * (r_vec/r) * (1 + 3*l_merc*l_merc/(r*r*c_sq)) / (r*r);
  }
  else{
    Fg = otherPlanet.mass * (r_vec/r) / (r*r);
  }
  return Fg;
}


vec Planet::distanceOther(Planet& otherPlanet, int i){
  //vec dr = pos - otherPlanet.pos;
  //return sqrt(dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2));

  //compute distance to other planet for all three directions at a given time (=i)
  vec dis(3);
  dis[0] = pos[i] - otherPlanet.pos[i];
  dis[1] = pos[i+N] - otherPlanet.pos[i+N];
  dis[2] = pos[i+2*N] - otherPlanet.pos[i+2*N];

  //return  sqrt(r_x*r_x + r_y*r_y + r_z*r_z);
  return dis;
}


vec Planet::gravitationalForce(Planet& otherPlanet, int i){

  vec r_vec =  - distanceOther(otherPlanet, i);
  double Beta = 2;
  double r = norm(r_vec);
  vec Fg;
  //if ((name == "Mercury" && otherPlanet.name = "Sun") || (name == "Sun" &&  otherPlanet.name = "Mercury")){
  if ( true == true){
    //cout << "Extra force in use" << endl;
    //cout << "orgForce" << otherPlanet.mass * (r_vec/r) / (pow(r,Beta)) << endl;
    //cout << "extraForce" << otherPlanet.mass * (r_vec/r) * (3*l_merc*l_merc/(r*r*c_sq)) / (pow(r,Beta)) << endl;
    Fg = otherPlanet.mass * (r_vec/r) * (1 + 3*l_merc*l_merc/(r*r*c_sq)) / (pow(r,Beta));
    //cout << "This planet: "<< name << ", f: " << Fg << endl;
  }
  else{
    Fg = otherPlanet.mass * (r_vec/r) / (pow(r,Beta));
  }
  return Fg;
}
