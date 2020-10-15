#include "solver.hpp"
using namespace std;
using namespace arma;

Solver::Solver(vector<Planet> sysPlanets, int N_val, double t_n_val, string sys){
  sysName = sys;
  planets = sysPlanets;
  N = N_val;
  t_n = t_n_val;
  //cout << "PRINT OUT PI:  " << pi << endl;
}

vec Solver::TotalAccelerationOnPlanet(Planet& planet, int index){ //calculate total acceleration on planet
  vec accel(3, fill::zeros); // acceleration vector [a_x, a_y, a_z] filled with zeros [0,0,0]
  for(int i =0; i < planets.size(); i++){ // find neighbour planets and calculate acceleration
    if(planets[i].name != planet.name){
      accel += planet.gravitationalForce(planets[i], index); //fetch gravitationalForce on planet due to neighbour planet for a given time (=index)

    }
  }
  return accel*G_scale;
}

void Solver::Euler(){
  method = "E";
  double h = t_n/N; //stepsize
    for (int j = 0; j < N-1; j++){
        for(int k=0; k < planets.size(); k++){ //for every planet compute position and velocity at specific time j
          vec h_accel = h*TotalAccelerationOnPlanet(planets[k], j); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
          planets[k].vel[j+1] =  planets[k].vel[j] + h_accel[0]; // update x velocity
          planets[k].vel[j+1+N] =  planets[k].vel[j+N] + h_accel[1]; // update y velocity
          planets[k].vel[j+1+2*N] =  planets[k].vel[j+2*N] + h_accel[2]; // update z velocity

          planets[k].pos[j+1] = planets[k].pos[j] + h*planets[k].vel[j]; // update x position
          planets[k].pos[j+1+N] = planets[k].pos[j+N] + h*planets[k].vel[j+N]; // update y position
          planets[k].pos[j+1+2*N] = planets[k].pos[j+2*N] + h*planets[k].vel[j+2*N]; // update z position
    }
  }
  return;

}


void Solver::EulerCromer(){
  method = "EC";
  double h = t_n/N; //stepsize
    for (int j = 0; j < N-1; j++){
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
  method = "VV";
  double h = t_n/N; //stepsize
  double h_2 = h/2.0; //stepsize
    for (int j = 0; j < N-1; j++){
        for(int k=0; k < planets.size(); k++){ //for every planet compute position and velocity at specific time j
          vec accel = TotalAccelerationOnPlanet(planets[k], j); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
          planets[k].pos[j+1] = planets[k].pos[j] + h*planets[k].vel[j] + h*h_2*accel[0]; // update x position
          planets[k].pos[j+1+N] = planets[k].pos[j+N] + h*planets[k].vel[j+N]+ h*h_2*accel[1]; // update y position
          planets[k].pos[j+1+2*N] = planets[k].pos[j+2*N] + h*planets[k].vel[j+2*N]+ h*h_2*accel[2]; // update z position
        }
        for(int k=0; k < planets.size(); k++){
          vec accel = TotalAccelerationOnPlanet(planets[k], j); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
          vec accel_next = TotalAccelerationOnPlanet(planets[k], j+1); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j+1
          planets[k].vel[j+1] =  planets[k].vel[j] + h_2*(accel[0]+accel_next[0]); // update x velocity
          planets[k].vel[j+1+N] =  planets[k].vel[j+N] + h_2*(accel[1]+accel_next[1]); // update y velocity
          planets[k].vel[j+1+2*N] =  planets[k].vel[j+2*N] + h_2*(accel[2]+accel_next[2]); // update z velocity

          //optimer med accel_next = accel
    }
  }

  return;
}

void Solver::testTotE(){
  double eps = 1e-7;
  //calculate total energy last and first timestep, check if diff smaller than eps
  double startE = 0;
  double endE = 0;
  for(int k = 0; k < planets.size() ; k++){
    startE += calcPE(k, 0) + calcKE(k, 0);
    endE += calcPE(k, N-1) + calcKE(k, N-1);
  }
  double error = abs(startE) - abs(endE);
  cout << "Start vs end energy: "<< startE << ", " << endE << endl;
  cout << "Error  = " << error << endl;
}

double Solver::calcKE(int k, int j){
  double K = 0;
  Planet plan = planets[k];
  double vx2 = pow(plan.vel[j], 2);
  double vy2 = pow(plan.vel[j+N], 2);
  double vz2 = pow(plan.vel[j+2*N], 2);
  K += 0.5 * plan.mass *(vx2 + vy2 + vz2);
  return K;
}

double Solver::calcPE(int k, int j){
  // Calculates potential energy for ONE planet by ONE timepoint
  double U = 0; //pot energy
  double m1 = planets[k].mass; //this mass
  double r; // distance between the 2 planets
  double m2; //other planet mass

  //Sum up potential energy for all the planets
  for(int i = 0; i < planets.size() ;i++)
    if( i != k){
      m2 = planets[i].mass;
      r = norm(planets[k].distanceOther(planets[i], j));
      U -= m1*m2*G_scale/r;
    }
  return U;
}

void Solver::VelocityVerlet2(){
  method = "VV2";
  double h = t_n/N;
  double h_2 = h/2;


}

void Solver::WriteResults(){
  ofstream ofile;
  ofile.open("Results/" + sysName + "_" + method + ".csv");
  ofile << method << ", " <<  t_n/N << ", " << t_n << endl;
  for (int j = 0; j < N; j ++){
    for (int k = 0; k < planets.size(); k ++){
      ofile <<
      planets[k].pos[j] << ", " <<
      planets[k].pos[j+N] << ", " <<
      planets[k].pos[j+2*N] << ", " <<
      planets[k].vel[j] << ", " <<
      planets[k].vel[j+N] << ", " <<
      planets[k].vel[j+2*N] << ", ";
    }
    ofile << endl;
  }
}

// double Solver::getTotalEnergy(){
//   //Gets energy of the whole system
//   double totE = 0;
//   for(int k=0; k < planets.size(); k++){
//     totE += planets[k].getPE();
//     totE += planets[k].getKE();
//   }
//   return totE;
// }
