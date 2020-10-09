#include "solver.hpp"
using namespace std;
using namespace arma;

Solver::Solver(vector<Planet> sysPlanets, int N_val, double t_n_val){
  planets = sysPlanets;
  N = N_val;
  t_n = t_n_val;
}

vec Solver::TotalAccelerationOnPlanet(Planet& planet, int index){ //calculate total acceleration on planet
  vec accel(3, fill::zeros); // acceleration vector [a_x, a_y, a_z] filled with zeros [0,0,0]
  for(int i =0; i < planets.size(); i++){ // find neighbour planets and calculate acceleration
    if(planets[i].name != planet.name){
      accel += planet.gravitationalForce(planets[i], index, N); //fetch gravitationalForce on planet due to neighbour planet for a given time (=index)
    } else{
      cout << "Calculate forces on " << planets[i].name << endl;
    }
  }
  return accel*G_scale;
}

void Solver::EulerCromer(){
  double h = t_n/N; //stepsize
    for (int j = 0; j < N; j++){
        for(int k=0; k < planets.size(); k++){ //for every planet compute position and velocity at specific time j
          vec h_accel = h*TotalAccelerationOnPlanet(planets[k], j); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
          planets[k].vel[j+1] =  planets[k].vel[j] + h_accel[0]; // update x velocity
          planets[k].vel[j+1+N] =  planets[k].vel[j+N] + h_accel[1]; // update y velocity
          planets[k].vel[j+1+2*N] =  planets[k].vel[j+2*N] + h_accel[2]; // update z velocity

          planets[k].pos[j+1] = planets[k].pos[j] + h*planets[k].vel[j+1]; // update x position
          planets[k].pos[j+1+N] = planets[k].pos[j+N] + h*planets[k].vel[j+1+N]; // update y position
          planets[k].pos[j+1+2*N] = planets[k].pos[j+2*N] + h*planets[k].vel[j+1+2*N]; // update z position
    }
  }
  return;
}


void Solver::VelocityVerlet(){
  double h = t_n/N; //stepsize
  double h_2 = h/2.0; //stepsize
    for (int j = 0; j < N; j++){
        for(int k=0; k < planets.size(); k++){ //for every planet compute position and velocity at specific time j
          vec accel = TotalAccelerationOnPlanet(planets[k], j); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
          //planets[k].pos[j+1] = planets[k].pos[j] + h*planets[k].vel[j] + h*h_2*accel[0]/2.0; // update x position
          //planets[k].pos[j+1+N] = planets[k].pos[j+N] + h*planets[k].vel[j+N]+ h*h_2*accel[1]/2.0; // update y position
          //planets[k].pos[j+1+2*N] = planets[k].pos[j+2*N] + h*planets[k].vel[j+2*N]+ h*h_2*accel[2]/2.0; // update z position

          //vec accel_next = TotalAccelerationOnPlanet(planets[k], j+1); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j+1
          //planets[k].vel[j+1] =  planets[k].vel[j] + h_2*(accel[0]+accel_next[0]); // update x velocity
          //planets[k].vel[j+1+N] =  planets[k].vel[j+N] + h_2*(accel[1]+accel_next[1]); // update y velocity
          //planets[k].vel[j+1+2*N] =  planets[k].vel[j+2*N] + h_2*(accel[2]+accel_next[2]); // update z velocity

          //optimer med accel_next = accel


    }
  }
  return;
}
