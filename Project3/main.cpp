#include "solver.hpp"
#include "planet.hpp"
#include <map>
using namespace arma;
using namespace std;

vector<Planet> read_initial(vector<string> object_names, int N_points);

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
  systems["systemC"] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
  systems["systemD"] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter"};

  //planets.push_back(Planet(1898.13e24/1988500e24, 5, 0., 0., 0., 2, 0., "Jupiter", N_points));
  vector<Planet> planets;
  planets = read_initial(systems[system], N_points);

  //planets.push_back(Planet(1, 0., 0., 0., 0., 0., 0., "Sun", N_points));
  //planets.push_back(Planet(5.97219e24/1988500e24, 1., 0., 0., 0., 1.42*2*pi, 0., "Earth", N_points));
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
  solv.WriteResults();
  //solv.VelocityVerlet();
  solv.testTotE();

  return 0;
}

vector<Planet> read_initial(vector<string> sys_names, int N_points){
  int N_objects = 10;
  vector<string> object_names = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
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
