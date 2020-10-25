#include "solver.hpp"
#include "planet.hpp"
#include <map>
using namespace arma;
using namespace std;

vector<Planet> read_initial(vector<string> object_names, int N_points);
vector<Planet> adjustedOrigin(vector<Planet> planets, int N);
int main(int argc, char const *argv[]) {
  string system = argv[1];
  string method = argv[2];
  int t_end = atoi(argv[3]);
  double h = stod(argv[4]);
  double pi = 2*acos(0.0);
  int N_points = (int) ( (double)t_end/h);
  //cout << N_points << endl;

  map<string, vector<string> > systems;
  systems["systemA"] = {"Sun", "Earth"};
  systems["systemB"] = {"Sun", "Earth", "Jupiter"};
  systems["systemC"] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
  systems["systemD"] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter"};
  systems["systemE"] = {"Sun", "Mercury"};

  //planets.push_back(Planet(1898.13e24/1988500e24, 5, 0., 0., 0., 2, 0., "Jupiter", N_points));
  vector<Planet> planets;
  //planets = read_initial(systems[system], N_points);
  planets.push_back(Planet(1, 0., 0., 0., 0., 0., 0., "Sun", N_points, 0));
  //planets.push_back(Planet(5.97219e24/1988500e24, 1., 0., 0., 0., 1.42*2*pi, 0., "Earth", N_points)); //testing excape velcotiy sqrt(2)*v_circular
  //planets.push_back(Planet(5.97219e24/1988500e24, 1, 0., 0., 0., 5, 0., "Earth", N_points)); //v_circular
  //planets.push_back(Planet(5.97219e24/1988500e24, 1., 0., 0., 0., 5, 0., "Earth", N_points)); //v_elliptical
  //planets.push_back(Planet(1.89813e27/1988500e24, 5.1, 0., 0., 0., 2*pi/sqrt(5.1), 0., "Jupiter", N_points)); //v_circular
  planets.push_back(Planet(3.285E23/1988500e24 , 0.3075, 0., 0., 0., 12.44, 0., "Mercury", N_points, 0));
  //planets.push_back(Planet(5.97219e24/1988500e24, 1., 0., 0., 0., 1.42*2*pi, 0., "Earth", N_points));
  planets = adjustedOrigin(planets, N_points);

  Solver solv(planets, N_points , t_end, system);

  if(method=="E"){
      solv.Euler();
  }
  if(method=="EC"){
      solv.EulerCromer();
  }
  if(method=="VV"){
      solv.VelocityVerlet();
  }

  if(method=="VV2"){
      solv.VertleNoStorage();
  }

  //solv.WriteResults();
  //solv.VelocityVerlet();
  //solv.testTotE();
  //solv.testAngMom();
  //solv.WritePeriResults();

  return 0;
}

vector<Planet> adjustedOrigin(vector<Planet> planets, int N){
  vec CoM(3);
  double Mtot;
  vec Moms(3); //Sum r * m
  vec MomTot(3); //Total momentum of system

  for (int i = 0; i < planets.size(); i++){
    Mtot += planets[i].mass;
    Moms(0) += planets[i].pos[0] * planets[i].mass;
    Moms(1) += planets[i].pos[N] * planets[i].mass;
    Moms(2) += planets[i].pos[2*N] * planets[i].mass;

    MomTot(0) += planets[i].vel[0] * planets[i].mass;
    MomTot(1) += planets[i].vel[N] * planets[i].mass;
    MomTot(2) += planets[i].vel[2*N] * planets[i].mass;
  }

  CoM = Moms/Mtot;


  for (int i = 0; i < planets.size(); i++){
    planets[i].pos[0] -= CoM(0);
    planets[i].pos[N] -= CoM(1);
    planets[i].pos[2*N] -= CoM(2);

    planets[i].vel[0] -= MomTot(0)/Mtot;
    planets[i].vel[N] -= MomTot(1)/Mtot;
    planets[i].vel[2*N] -= MomTot(2)/Mtot;

  }
  return planets;
}



vector<Planet> read_initial(vector<string> sys_names, int N_points){
  int N_objects = 10;
  vector<string> object_names = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
  double *x, *y, *z, *vx, *vy, *vz; //To store initial conditions for each particle.
  double *mass; //Store mass of particles.
  x = new double[N_objects];
  y = new double[N_objects];
  z = new double[N_objects];
  vx = new double[N_objects];
  vy = new double[N_objects];
  vz = new double[N_objects];
  mass = new double[N_objects];

  const char* filename_pos_and_vel = "pos_vel_initial.txt";   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  const char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.

  //Open files
  FILE *fp_init = fopen(filename_pos_and_vel, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r"); //Open file to read.

  //Loop over each particle and extract its mass and initial conditions:
  for (int i = 0; i < N_objects; i++){
  	fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]); // One %lf (lf=long float or double) for each floating point number on each line of the file.
  	fscanf(fp_mass, "%lf", &mass[i]); //Extract mass for particle i.
  }

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.



  vector<Planet> planets;
  for (int i = 0; i < object_names.size(); i ++){
    for (int j = 0; j < sys_names.size(); j ++){
      if (object_names[i] == sys_names[j])
        planets.push_back(Planet(mass[i]/mass[0], x[i], y[i], z[i], vx[i]*365, vy[i]*365, vz[i]*365, sys_names[j], N_points));
    }
  }

  return  planets;
}
