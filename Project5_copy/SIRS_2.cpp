#include "SIRS.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <armadillo>
using namespace std;
using namespace arma;

SIRS::SIRS(double S_, double I_, double a_, double b_, double c_, double t_){
  y = vec(3);
  y(0) = S_;
  y(1) = I_;
  y(2) = 0;

  S0 = S_;
  I0 = I_;
  R0 = 0;

  N = S_+ I_;
  a = a_;
  b = b_;
  c = c_;

  t = t_;
}
void SIRS::specRK4(double dt_){
  dt = dt_;
  dy = vec(3);
  useDer1 = false;
  num_pts = int(t/dt);
}

void SIRS::specMC(int MC_cyc){ //MC
  MC_cycles = MC_cyc;
  useVD = false;
}

void SIRS::specRK4_VD(double dt_, double e_, double d_, double dI){ //RK4 with e and d
 dt = dt_;
 dy = vec(3);
 d = d_;
 e = e_;
 d_I = dI;
 num_pts = int(t/dt);
 useDer1 = false;
}

void SIRS::specMC_VD(int MC_cyc, double e_, double d_, double dI){ //MC with e and d
  MC_cycles = MC_cyc;
  d = d_;
  e = e_;
  d_I = dI;
  useVD = true;
}


vec SIRS::derivatives(vec yt){
  // S, I, R = y(0), y(1), y(2)
  dy(0) = c*yt(2) - a*yt(1)*yt(0)/N;
  dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  return dy;
}

vec SIRS::derivatives2(vec yt){
  //Update N since people die and get born
  dy(0) = c*yt(2) - yt(0)*(a*yt(1)/N+d)+e*N;
  dy(1) = yt(1)*(a*yt(0)/N - b-d-d_I);
  dy(2) = b*yt(1) - (c+d)*yt(2);
  return dy;
}


void SIRS::solveRK4(string filename){
  ofstream ofile;
  ofile.open(filename + "_RK4.csv");

  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", RK4" << ", std" <<  endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
  //for (double i = 0; i < t; i += dt){
  cout << "Num pts in RK4 solve: " <<  num_pts << endl;
  for(int i = 0; i < num_pts; i++){
    //cout << i << endl;
    rk4(useDer1);
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

  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", MC"  << ", std" << endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
  cout << "dt: " << dt << endl;
  //Kjøre MC, ca 1000 ganger. Ha 3 t/dt lang array med verdier for S,I,R.
  //Ta gjennomsnitt av alle kjøringene og skriv dette til fil.
  //Må vel resetee S0 etc hver gang vel, starte fra t0 for hver eksperiment,
  //reset initial conditions with func
  bool useVD = false;
  double progress = 0;
  for (int j = 0; j < MC_cycles; j ++){
    if(j+1 >= progress){
      cout << progress*100/MC_cycles << "% done." << endl; //Interval time: " << endl;
      progress += 0.1*MC_cycles;
    }
      //stop = clock();
      //double timeInterval = ( (stop - start)/(double)CLOCKS_PER_SEC );
      //totTime += timeInterval;
      //start = clock();
    for (int i = 0; i < num_pts; i++){
      /*
      if(i == 2){
        cout << "MC = " << j << endl;
        cout << S0 << ", " << ", " << I0 << ", " << R0 << endl;
        cout << y(0) << endl;
        cout << S_mc(0) << endl;
        cout << S_mc(1) << "\n" <<  endl;
      }*/
      MonteCarlo();
      S_mc(i) += y(0);
      I_mc(i) += y(1);
      R_mc(i) += y(2);
    }
    reset_states();
  }
  S_mc = S_mc/MC_cycles;
  I_mc = I_mc/MC_cycles;
  R_mc = R_mc/MC_cycles;

  for (int i = 0; i < num_pts; i++){
    //cout << S_mc(i) << ", " << I_mc(i) << ", " << R_mc(i) << endl;
    ofile << S_mc(i) << ", " << I_mc(i) << ", " << R_mc(i) << endl;
  }
}


void SIRS::rk4(bool useDer1){
  vec K1(3), K2(3), K3(3), K4(3);

  if(useDer1 == true){
    //cout << "Using derivate1 function." << endl;
    K1 = dt*derivatives(y);
    //cout << "K1: " << K1 << endl;
    K2 = dt*derivatives(y+0.5*K1);
    K3 = dt*derivatives(y+0.5*K2);
    K4 = dt*derivatives(y+K3);
  }

  else if(useDer1 == false){
    //cout << "Using der2" << endl;
    K1 = dt*derivatives2(y);
    //cout << "K1: " << K1 << endl;
    K2 = dt*derivatives2(y+0.5*K1);
    K3 = dt*derivatives2(y+0.5*K2);
    K4 = dt*derivatives2(y+K3);
  }
  //cout << K2 << endl;
  y = y + (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0;
  y(2) = N - y(1) - y(0);
  //cout << "RK4" << endl;
  //cout <<"S: " << y(0) << ", I: "<< y(1) << endl;
}


void SIRS::MonteCarlo(){
  // S, I, R = y(0), y(1), y(2)
  //Tansition probabilities and dt
  pR_S = (c*y(2))*dt;
  pS_I = (a*y(0)*y(1)/N)*dt;
  pI_R = (b*y(1))*dt;

  int RS_count = 0, SI_count = 0, IR_count = 0;
  int bornS = 0, diedS = 0, diedI = 0, diedI_disease = 0, diedR = 0;

  r = rand() % 100001;
  if (r/100000 < pR_S)
      RS_count ++;

  r = rand() % 100001;
  if (r/100000 < pS_I)
      SI_count ++;

  r = rand() % 100001;
  if (r/100000 < pI_R)
      IR_count ++;

  if(useVD == true){
    r = rand() % 100001;
    if (r/100000 < e*N*dt) //Is birth rate given for dt = 1??
        bornS ++;

    r = rand() % 100001;
    if (r/100000 < d*y(0)*dt) //Is birth rate given for dt = 1??
        diedS ++;

    r = rand() % 100001;
    if (r/100000 < d*y(1)*dt) //Is birth rate given for dt = 1??
        diedI ++;

    r = rand() % 100001;
    if (r/100000 < d*y(1)*dt) //Is birth rate given for dt = 1??
        diedI_disease ++;

    r = rand() % 100001;
    if (r/100000 < d*y(2)*dt) //Is birth rate given for dt = 1??
        diedR ++;
    N = y(0) + y(1) + y(2); // Update tot pop after deaths and births
    // Do we have to update new dt since N change?
  }
  y(0) += RS_count - SI_count + bornS - diedS;
  y(1) += SI_count - IR_count - diedI - diedI_disease;
  y(2) += IR_count - RS_count - diedR;
  //cout << "MC" << endl;
  //cout <<"S: " << y(0) << ", I: "<< y(1) << endl;
}




// inline void writeResults(fstream file);
//   file << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
