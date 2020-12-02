#include "SIRS.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;


SIRS::SIRS(double S_, double I_, double a_, double b_, double c_, double t_, double dt_){
  y = vec(3);
  dy = vec(3);
  y(0)=S_;
  y(1) = I_;
  y(2) = 0;

  N = S_+ I_;
  a = a_;
  b = b_;
  c = c_;

  t = t_;
  dt = dt_;
  dt2 = dt/2;
}

vec SIRS::derivatives(vec yt){
  dy(0) = c*yt(2) - a*yt(1)*yt(0)/N;
  dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  return dy;
}


void SIRS::solve(string filename, string method){
  ofstream ofile;
  ofile.open(filename + "_" + method +".csv");

  if(method == "RK4"){
    ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", " << method << endl;
    ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
    for (double i = 0; i < t; i += dt){
      rk4();
      ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
      //writeResults(ofile);
    }
  }
  else if(method == "MC"){
    dt = min({4/a/N , 1/b/N, 1/c/N});
    ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", " << method << endl;
    ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
    cout << "dt: " << dt << endl;
    for (double i = 0; i < t; i += dt){
      MonteCarlo();
      ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
      //writeResults(ofile);
    }
  }

  cout << "Ended solve" << endl;
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
  double pR_S  = (c*y(2))*dt;
  double pS_I = (a*y(0)*y(1)/N)*dt;
  double pI_R = (b*y(1))*dt;
  double r;

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
