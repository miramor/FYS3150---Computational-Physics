#include "solver.hpp"
#include "time.h"
using namespace std;
using namespace arma;

Solver::Solver(vector<Planet> sysPlanets, int N_val, double t_n_val, string sys){
  sysName = sys;
  planets = sysPlanets;
  N = N_val;
  t_n = t_n_val;

  cout << "\n" << sys << " created, contains:" << endl;
  for(int k = 0; k < planets.size(); k ++){
    cout << planets[k].name << endl;
  }
  cout << "\n" << endl;

  //Pre calculate ang momentum for Mercury, but need to find a way to use it in Planet
  if(sys == "systemE"){
    Planet sun = planets[0];
    Planet merc = planets[1];
    vec r_vec =  - planets[1].distanceOther_opt(planets[0], 0);
    vec v_vec(3);

    v_vec[0] = merc.vel[0]-sun.vel[0];
    v_vec[1] = merc.vel[1]-sun.vel[1];
    v_vec[2] = merc.vel[2]-sun.vel[2];

    double l_merc = norm(cross(r_vec,v_vec)); //Angular orbital momentum, only calculate once
    //cout << "ang momemnt merc:  " << l_merc << endl;
    planets[1].l_merc = l_merc; //give Mercury access to this to be used for the additional force
  }
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

  double progress = 0.1;
  clock_t start, stop;
  double totTime = 0;
  start = clock();
    for (int j = 0; j < N-1; j++){

      if(j >= progress*(N-2)){
        stop = clock();
        double timeInterval = ( (stop - start)/(double)CLOCKS_PER_SEC );
        totTime += timeInterval;
        cout << progress*100 << "% done. Interval time: " << timeInterval << endl;
        progress = progress + 0.1;
        start = clock();
      }

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
    }
  }

  cout << "Total time(s)--time(m): " << totTime << " -- " << totTime/60 << endl;
  return;
}


void Solver::testAngMom(){
  double eps = 1e-5;

  double l1 = calcL(0);
  double l2 = calcL(N-1);

  double error = l1 - l2;

  //cout << "Error in Angular Momenutm: " << error << endl;
  cout << "Relative error AngMom: " << error/l1 << endl;
}

double Solver::calcL(int index){
  vec r_vec = planets[0].distanceOther(planets[1], index);

  vec v_vec(3);
  v_vec[0] = planets[0].vel[index]-planets[1].vel[index];
  v_vec[1] = planets[0].vel[index+N]-planets[1].vel[index+N];
  v_vec[2] = planets[0].vel[index+2*N]-planets[1].vel[index+2*N];

  double l = norm(cross(r_vec,v_vec));
  return l;

}


void Solver::testTotE(){
  double eps = 1e-5;
  //calculate total energy last and first timestep, check if diff smaller than eps
  double startE = 0;
  double endE = 0;
  for(int k = 0; k < planets.size() ; k++){
    startE += calcPE(k, 0)/2 + calcKE(k, 0);
    endE += calcPE(k, N-1)/2 + calcKE(k, N-1);
  }
  double error = abs(startE) - abs(endE);
  /*
  if (abs(error) > eps)
      cout << "The total energy is not conserved.\n Error: " << abs(error) << endl;
  */
  //cout << "Error in total energy: " << error << endl;
  cout << "Relative change in total energy: " << error/startE << endl;
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
      //Use if central gravitational force is changed to not being inverse-square
      // double beta = 2;
      // U += m1*m2*G_scale/((beta-1)*pow(r,beta-1);
    }
  return U;
}


void Solver::WriteResults(){
  ofstream ofile;
  ofile.open("Results/" + sysName + "_" + method + ".csv");
  ofile << method << ", " <<  t_n/N << ", " << t_n << endl;
  for (int j = 0; j < N; j += 10){
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

vec Solver::TotalAccelerationOnPlanet_opt(Planet& planet, bool useCurr){ //calculate total acceleration on planet
  vec accel(3, fill::zeros); // acceleration vector [a_x, a_y, a_z] filled with zeros [0,0,0]

  for(int i =0; i < planets.size(); i++){ // find neighbour planets and calculate acceleration
    if(planets[i].name != planet.name){
      accel += planet.gravitationalForce_opt(planets[i], useCurr); //fetch gravitationalForce on planet due to neighbour planet for a given time (=index)

    }
  }
  return accel*G_scale;
}

void Solver::VertleNoStorage(){
  /*
    Metod for solving VelocityVerlet, but without storing all the data inside the planet
    object and a full csv file. Writes the relevant data to "Results/Peri_Results.csv".
    Format [x, y, r, index(j)].

    Only used for systemE, containing Mercury & Sun. Uses alternate methods in class
    indicated by the name of the method ending with "_opt".
  */
  method = "VV2";
  double h = t_n/N; //stepsize
  double h_2 = h/2.0; //stepsize
  //ofstream ofile;
  //ofile.open("Results/" + sysName + "_" + method + ".csv");
  //ofile << "relevant data for perihelion, mercury & sun:" << endl;

  /*
  for(int k=0; k < planets.size(); k++){
    Planet plan = planets[k];
    ofile << plan.pos[0] << " ,"  << plan.pos[1] << " ,"  << plan.pos[2] << ", ";
    ofile << plan.vel[0] << " ,"  << plan.vel[1] << " ,"  << plan.vel[2] << ", ";
  }*/

  //ofile << endl;
  double progress = 0.1;
  clock_t start, stop;
  double totTime = 0;
  start = clock();

  ofstream ofilePeri;
  ofilePeri.open("Results/Peri_Results.csv");

  double r0,r1,r2;
  vec r0_vec(3), r1_vec(3), r2_vec(3);
  r0_vec = planets[1].distanceOther_opt(planets[0], false);

  ofilePeri << setprecision(20);
  ofilePeri << r0_vec[0] << " , "<< r0_vec[1] << " , " << r0_vec[2] << " , " << norm(r0_vec) << " , " << 0 << endl;

    //Start calculations for all timesteps
    for (int j = 0; j < N-1; j++){

      // Prints out progress
      if(j >= progress*(N-2)){
        stop = clock();
        double timeInterval = ( (stop - start)/(double)CLOCKS_PER_SEC );
        totTime += timeInterval;
        cout << progress*100 << "% done. Interval time: " << timeInterval << endl;
        progress = progress + 0.1;
        start = clock();
      }

      for(int k=0; k < planets.size(); k++){ //for every planet compute position and velocity at specific time j
        vec accel = TotalAccelerationOnPlanet_opt(planets[k], false); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
        planets[k].pos[3] = planets[k].pos[0] + h*planets[k].vel[0] + h*h_2*accel[0]; // update x position
        planets[k].pos[4] = planets[k].pos[1] + h*planets[k].vel[1]+ h*h_2*accel[1]; // update y position
        planets[k].pos[5] = planets[k].pos[2] + h*planets[k].vel[2]+ h*h_2*accel[2]; // update z position
      }


      for(int k=0; k < planets.size(); k++){
        vec accel = TotalAccelerationOnPlanet_opt(planets[k], false); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j
        vec accel_next = TotalAccelerationOnPlanet_opt(planets[k], true); //fetch planet's acceleration vector [a_x, a_y, a_z] times h at a given time j+1

        planets[k].vel[3] =  planets[k].vel[0] + h_2*(accel[0]+accel_next[0]); // update x velocity
        planets[k].vel[4] =  planets[k].vel[1] + h_2*(accel[1]+accel_next[1]); // update y velocity
        planets[k].vel[5] =  planets[k].vel[2] + h_2*(accel[2]+accel_next[2]); // update z velocity

      }

      if (j == 0){
        r1_vec = planets[1].distanceOther_opt(planets[0], true);
        r1 = norm(r1_vec);

      }
      else{
        r2_vec = planets[1].distanceOther_opt(planets[0],true);
        r2 = norm(r2_vec);

        //Finds the point where planet is closest to the sun and writes to file.
        if (r0 > r1 && r1 < r2){
          ofilePeri << r1_vec[0] << " , "<< r1_vec[1] << " , " << r1_vec[2] << r1 << " , " << j << endl;
        }
        //Update for new timestep
        r0_vec = r1_vec;
        r0 = norm(r0_vec);
        r1_vec = r2_vec;
        r1 = norm(r1_vec);
      }


      //Write out new results
      /*
      if(j % 1000){
        for(int k=0; k < planets.size(); k++){
          Planet plan = planets[k];
          ofile << plan.pos[3] << " ,"  << plan.pos[4] << " ,"  << plan.pos[5] << ", ";
          ofile << plan.vel[3] << " ,"  << plan.vel[4] << " ,"  << plan.vel[5] << ", ";

        }
      ofile << endl;
      }
      */

      for(int k=0; k < planets.size(); k++){
        planets[k].pos[0] = planets[k].pos[3];
        planets[k].pos[1] = planets[k].pos[4];
        planets[k].pos[2] = planets[k].pos[5];

        planets[k].vel[0] = planets[k].vel[3];
        planets[k].vel[1] = planets[k].vel[4];
        planets[k].vel[2] = planets[k].vel[5];
      }


    }
  return;
}
