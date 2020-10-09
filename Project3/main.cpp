#include "solver.hpp"
#include "planet.hpp"

using namespace arma;
using namespace std;

vector<Planet> read_initial(int N_objects, int N_points, vector<string> object_names);

int main(int argc, char const *argv[]) {

  int N_points = 10000;
  Planet sun(1, 0., 0., 0., 0., 0., 0., "sun", N_points);
  Planet earth(5.972e24/1.989e30, 1., 0., 0., 0., 2*3.141592, 0., "earth", N_points);

  vector<Planet> planets;
//  planets = read_initial(10, N_points)
  planets.push_back(sun);
  planets.push_back(earth);

  Solver solv(planets, N_points , 1, "systemA");
  solv.EulerCromer();


  //Planet jupyter =planet(1.898E27 kg, vec p(3) = {2, 0, 0}, vec v(3) = {0, 1, 0});

  //Vector<Planet> planets;
  //Vector<string> object_names = {"Sun", "Earth"};
  //planets = read_initial(10, 100, object_names);
  return 0;
}
/*
vector<Planet> read_initial(int N_objects, int N_points, vec<string> object_names){
  double *x, *y, *z, *vx, *vy, *vz; //To store initial conditions for each particle.
  double *mass; //Store mass of particles.
  string *names; //Store names of objects.
  x = new double[N_objects];
  y = new double[N_objects];
  z = new double[N_objects];
  vx = new double[N_objects];
  vy = new double[N_objects];
  vz = new double[N_objects];
  mass = new double[N_objects];
  names = new string[N_objects];

  char* filename_pos_and_vel = "pos_vel_initial.txt";   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.

  //Open files
  FILE *fp_init = fopen(filename_pos_and_vel, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r"); //Open file to read.

  //Loop over each particle and extract its mass and initial conditions:
  for (int i = 0; i < N_objects; i++){
  	fscanf(fp_init, "%s %lf %lf %lf %lf %lf %lf", &names[i], &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]); // One %lf (lf=long float or double) for each floating point number on each line of the file.
  	fscanf(fp_mass, "%s %lf", &names[i], &mass[i]); //Extract mass for particle i.
  }

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.

  vec<Planet> planets[object_names.size()];
  for (int i = 0; i < N_objects; i++){
    for (int j = 0; j < object_names.size(); j++){
      if (names[i] == object_names[j]){
        planets[i] = Planet(names[i], mass[i], {x[i], y[i], z[i]}, {vx[i], vy[i], vz[i]});
      }
    }
  }

  return planets;
}
*/
