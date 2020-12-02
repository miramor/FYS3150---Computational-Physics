#include "SIRS.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;


SIRS::SIRS(double S_, double I_, double a_, double b_, double c_, double t_, double dt_){

  // y.push_back(S_);
  // y.push_back(I_);
  // y.push_back(0);

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
  // dy(0) = c*(N-yt(0)-yt(1)) - a*yt(1)*yt(0)/N;
  // dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  //cout<<"dy: " << dy << endl;
  return dy;
}


void SIRS::solve(string filename){

  ofstream ofile;
  ofile.open(filename + ".csv");
  //vector <double> y, dy, ynext;
  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
  //t = 2;


  for (double i = 0; i < t; i += dt){
    rk4();
    ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
    //writeResults(ofile);
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

// inline void writeResults(fstream file);
//   file << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
