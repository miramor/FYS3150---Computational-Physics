#include "SIRS.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;

SIRS::SIRS(double S_, double I_, double a_, double b_, double c_, double t_, double dt_){ //RK4
  y = vec(3);
  dy = vec(3);
  y(0) = S_;
  y(1) = I_;
  y(2) = 0;

  N = S_+ I_;
  a = a_;
  b = b_;
  c = c_;

  t = t_;
  dt = dt_;
  num_pts = int(t/dt);
}

SIRS::SIRS(double S_, double I_, double a_, double b_, double c_, double t_, int MC_cyc){ //MC
  MC_cycles = MC_cyc;
  y = vec(3);

  S0 = S_;
  I0 = I_;
  R0 = 0;

  y(0) = S_;
  y(1) = I_;
  y(2) = 0;

  N = S_+ I_;
  a = a_;
  b = b_;
  c = c_;

  t = t_;
}


vec SIRS::derivatives(vec yt){
  dy(0) = c*yt(2) - a*yt(1)*yt(0)/N;
  dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  return dy;
}


void SIRS::solveRK4(string filename){
  ofstream ofile;
  ofile.open(filename + "_RK4.csv");

  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", RK4" << endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
  for (double i = 0; i < t; i += dt){
    rk4();
    ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
    //writeResults(ofile);
  }
}

void SIRS::reset_states(){
  y(0) = S0;
  y(1) = I0;
  y(2) = R0;
}

void SIRS::solveMC(string filename){
  ofstream ofile;
  ofile.open(filename + "_MC.csv");

  dt = min({4/a/N , 1/b/N, 1/c/N});
  num_pts = int(t/dt);

  S_mc = vec(num_pts, fill::zeros);
  I_mc = vec(num_pts, fill::zeros);
  R_mc = vec(num_pts, fill::zeros);

  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", MC"  << endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
  cout << "dt: " << dt << endl;
  //Kjøre MC, ca 1000 ganger. Ha 3 t/dt lang array med verdier for S,I,R.
  //Ta gjennomsnitt av alle kjøringene og skriv dette til fil.
  //Må vel resetee S0 etc hver gang vel, starte fra t0 for hver eksperiment,
  //reset initial conditions with func
  double progress = 0;
  for (int j = 0; j < MC_cycles; j ++){
    if( ){
      cout << progress << "% done. Interval time: " << endl;

    }
      //stop = clock();
      //double timeInterval = ( (stop - start)/(double)CLOCKS_PER_SEC );
      //totTime += timeInterval;
      //start = clock();
    }

    for (int i = 0; i < num_pts; i++){
      MonteCarlo();
      S_mc(i) += y(0);
      I_mc(i) += y(1);
      R_mc(i) += y(2);
    }
    reset_states();
    progress += 100/MC_cycles;
  }
  S_mc = S_mc/MC_cycles;
  I_mc = I_mc/MC_cycles;
  R_mc = R_mc/MC_cycles;

  for (int i = 0; i < num_pts; i++){
    ofile << S_mc(i) << ", " << I_mc(i) << ", " << R_mc(i) << endl;
  }
}


void SIRS::rk4(){
  vec K1(3), K2(3), K3(3), K4(3);

  K1 = dt*derivatives(y);
  //cout << "K1: " << K1 << endl;
  K2 = dt*derivatives(y+0.5*K1);
  K3 = dt*derivatives(y+0.5*K2);
  K4 = dt*derivatives(y+K3);
  //cout << K2 << endl;
  y = y + (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0;
  y(2) = N - y(1) - y(0);
}


void SIRS::MonteCarlo(){
  //Tansition probabilities and dt
  pR_S  = (c*y(2))*dt;
  pS_I = (a*y(0)*y(1)/N)*dt;
  pI_R = (b*y(1))*dt;

  int S_count = 0, I_count = 0, R_count = 0;


  r = rand() % 100001;
  if (r/100000 < pR_S)
      S_count ++;

  r = rand() % 100001;
  if (r/100000 < pS_I)
      I_count ++;

  r = rand() % 100001;
  if (r/100000 < pI_R)
      R_count ++;

  //cout << S_count << " " << I_count  << " " << R_count << endl;
  y(0) += R_count - S_count;
  y(1) += S_count - I_count;
  y(2) += I_count - R_count;

}




// inline void writeResults(fstream file);
//   file << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
