#include "solver.hpp"
using namespace std;
using namespace arma;


Solver::TotalForceOnPlanet(const Planet& planet){ //calculate total force on planet
  double force = 0;
  for(int i =0; i < planets.size(); i++){
    if(planets[i] != planet){
      force += planet.gravitationForce(planets[i]) //fetch gravitationalForce between the two planets
    }else{
      cout << "Calculate forces on " << planets[i] << endl;
    }
  }
  return force*G_scale*planet.mass
}

Solver::EulerCromer(int N){

}

// use for Euler-Cromer method 
for (int i = 0; i < M; i++){
    for (int j = 0; j < M; j++){
      for (int k = 0; k < N; k++){
        C[i*M + j] += A[i*N + k]*B[k*M + j];
      }
    }
  }
