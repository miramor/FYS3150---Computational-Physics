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


  //Contains all planets used when solving
  vector<Planet> planets;

  int dataChoice;
  cout << "Choose Custom or NASA data: (1)Custom, (2)NASA" << endl;
  cin >>  dataChoice;


  map<string, vector<string> > systems;
  systems["systemA"] = {"Sun", "Earth"};
  systems["systemB"] = {"Sun", "Earth", "Jupiter"};
  systems["systemC"] = {"Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
  systems["systemE"] = {"Sun", "Mercury"};

  if(system == "systemE" && method == "VV2" && dataChoice == 2){
    cout << "Cannot run VV2 for NASA data, changes to custom." << endl;
    dataChoice = 1;
  }

  if(dataChoice == 2){
    planets = read_initial(systems[system], N_points);
  }

  else{

    Planet sun = Planet(1, 0., 0., 0., 0., 0., 0., "Sun", N_points);

    if (system == "systemA") {
      vector<double> velocities = {2*pi, 5, 1.42*2*pi};
      int choice;
      cout << "Choose orbit for earth: circular(0), eliptical(1), escape vel(2)" << endl;
      cin >> choice;

      Planet earth = Planet(5.97219e24/1988500e24, 1, 0., 0., 0., velocities[choice], 0., "Earth", N_points);
      planets.push_back(sun);
      planets.push_back(earth);
    }

    if (system == "systemB") {
      Planet earth = Planet(5.97219e24/1988500e24, 1, 0., 0., 0., 2*pi, 0., "Earth", N_points);
      vector<double> massScale = {1, 10, 1000};
      int choice;
      cout << "Choose weigthed Jupiter: regular(0), times10(1), times1000(2)" << endl;
      cin >> choice;

      Planet jupiter = Planet(massScale[choice]*1.89813e27/1988500e24, 5.1, 0., 0., 0., 2*pi/sqrt(5.1), 0., "Jupiter", N_points);
      planets.push_back(sun);
      planets.push_back(earth); //earth circular
      planets.push_back(jupiter);//adjusted for list
    }

    if (system == "systemC") {
      cout << "Cant make custom for this system. Runs with NASA data instead" << endl;
      planets = read_initial(systems[system], N_points);
    }

    if (system == "systemE") {
      Planet sun = Planet(1, 0., 0., 0., 0., 0., 0., "Sun", N_points, 0);
      Planet mercury_opt = Planet(3.285E23/1988500e24 , 0.3075, 0., 0., 0., 12.44, 0., "Mercury", N_points, 0);
      planets.push_back(sun);
      planets.push_back(mercury_opt);

      cout << "Can only use VV2, runs with method VV2" << endl;
      method = "VV2";
    }
  }


  if( system != "systemE"){
    int choiceRef;
    cout << "\nEnter 1 for adjusted frame of reference to Center of Mass, any other number to not." << endl;
    cin >> choiceRef;

    if(choiceRef == 1){
      planets = adjustedOrigin(planets, N_points);
    }
  }

  Solver solv(planets, N_points , t_end, system);
  if(method=="E"){
    cout << "****************************" << endl;
    cout << "Euler started" << endl;
    solv.Euler();
    cout << "Finished Euler" << endl;
    cout << "****************************\n" << endl;
  }
  if(method=="EC"){
    cout << "****************************" << endl;
    cout << "EulerCromer started" << endl;
    solv.EulerCromer();
    cout << "Finished EulerCromer" << endl;
    cout << "****************************\n" << endl;
  }
  if(method=="VV"){
    cout << "****************************" << endl;
    cout << "VV started" << endl;
    solv.VelocityVerlet();
    cout << "Finished VV" << endl;
    cout << "****************************\n" << endl;
  }

  if(method=="VV2"){
    cout << "****************************" << endl;
    cout << "VV2 started" << endl;
    solv.VertleNoStorage();
    cout << "Finished VV2" << endl;
    cout << "****************************\n" << endl;
    return 0;
  }

  solv.WriteResults();

  solv.testTotE();
  solv.testAngMom();

  return 0;
}

vector<Planet> adjustedOrigin(vector<Planet> planets, int N){
  double Mtot = 0;
  vec CoM(3, fill::zeros), Vc(3, fill::zeros), Moms(3, fill::zeros), MomTot(3, fill::zeros);

  //Calculating total momentum of system
  for (int i = 0; i < planets.size(); i++){
    Mtot += planets[i].mass;
    Moms(0) += planets[i].pos[0] * planets[i].mass;
    Moms(1) += planets[i].pos[N] * planets[i].mass;
    Moms(2) += planets[i].pos[2*N] * planets[i].mass;

    MomTot(0) += planets[i].vel[0] * planets[i].mass;
    MomTot(1) += planets[i].vel[N] * planets[i].mass;
    MomTot(2) += planets[i].vel[2*N] * planets[i].mass;
  }

  CoM = Moms/Mtot; //Center of mass
  Vc = MomTot/Mtot; //Center of mass velocity


  //Adjusting position and velocity such that Center of Mass it at rest
  for (int i = 0; i < planets.size(); i++){
    planets[i].pos[0] -= CoM(0);
    planets[i].pos[N] -= CoM(1);
    planets[i].pos[2*N] -= CoM(2);

    planets[i].vel[0] -= Vc(0);
    planets[i].vel[N] -= Vc(1);
    planets[i].vel[2*N] -= Vc(2);

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
