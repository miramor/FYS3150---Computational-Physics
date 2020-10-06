#include "solver.hpp"
using namespace std;
using namespace arma;


Solver::TotalAccelerationOnPlanet(const Planet& planet){ //calculate total force on planet
  double accel = 0;
  for(int i =0; i < planets.size(); i++){
    if(planets[i] != planet){
      accel += planet.gravitationForce(planets[i]) //fetch gravitationalForce on planet due to neighbour planet
    } else{
      cout << "Calculate forces on " << planets[i] << endl;
    }
  }
  return accel*G_scale
}

Solver::EulerCromer(int N, double t_end, ve<Planet> planets){
  double h = t_end/N; //stepsize
    for (int j = 0; j < N; j++){
        for(k=0; k < planets.size(); k++){ //for every planet compute position and velocity at specific time
          double h_accel = h*TotalAccelerationOnPlanet(planet[k]) //planet's acceleration times h
          planets[k].vel[j+1] = planets.vel[j] + h_accel; // update x velocity
          planets[k].vel[j+1+N] = planets.vel[j+N] + h_accel; // update y velocity
          planets[k].vel[j+1+2*N] = planets.vel[j+2*N] + h_accel; // update z velocity

          planets[k].pos[j+1] = planets[k].pos[j] + h*planets[k].vel[j+1]; // update x position
          planets[k].pos[j+1+N] = planets[k].pos[j+N] + h*planets[k].vel[j+1+N]; // update y position
          planets[k].pos[j+1+2*N] = planets[k].pos[j+2*N] + h*planets[k].vel[j+1+2*N]; // update z position
    }
  }
}
