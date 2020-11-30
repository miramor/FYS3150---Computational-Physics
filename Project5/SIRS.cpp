#include "SIRS.hpp"
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

SIRS::SIRS(double S_, double I_, double a_, double b_, double c_, double t_, double dt_){

  y.push_back(S_);
  y.push_back(I_):
  y.push_back(0);

  dy.resize(3);
  N = S_+ I_;
  a = a_;
  b = b_;
  c = c_;

  t = t_;
  dt = dt_;
  dt2 = dt/2;
}

vector <double> SIRS::derivatives(yt){
  dy[0] = c*yt[2] - a*yt[1]*yt[0]/N;
  dy[1] = a*yt[0]*yt[0]/N - b*yt[1];
  return dy;
}


void SIRS::solve(string filename){
  fstream ofile;
  ofile.open(filename);
  vector <double> y, dy, ynext;

  for (int i = 0; i < t; i ++){
    rk4();
    writeResults(ofile);
  }

}

void SIRS::rk4(){
  vector <double> K1(3), K2(3), K3(3), K4(3);

  K1 = dt*derivatives(y);
  K2 = dt2*derivatives(y+K1);
  K3 = dt2*derivatives(y+K2);
  K4 = dt*derivatives(y+K3);
  cout << K2 << endl;
  y = y + (K1 + 2.*K2 + 2.*K3 + K4)/6;
  y[2] = N - y[1] - y[0];
}

inline void writeResults(fstream file);
  file << y[0] << ", " << y[1] << ", " <<  y[2] << endl;
