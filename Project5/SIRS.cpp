#include "SIRS.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <math.h>       /* cos */
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
  useVD = false;
  num_pts = int(t/dt);
}

void SIRS::specRK4(double dt_, double f_){
  dt = dt_;
  cout << useV << endl;
  dy = vec(3);
  useVD = false;
  num_pts = int(t/dt);
  useV = true;
  f = f_;

}

void SIRS::specMC(int MC_cyc){ //MC
  MC_cycles = MC_cyc;
  useVD = false;
}

void SIRS::specMC(int MC_cyc, double f_){ //MC
  MC_cycles = MC_cyc;
  useVD = false;
  useV = true;
  f = f_;
}

void SIRS::specRK4_VD(double dt_, double e_, double d_, double dI){ //RK4 with e and d
 dt = dt_;
 dy = vec(3);
 d = d_;
 e = e_;
 d_I = dI;
 num_pts = int(t/dt);
 useVD = true;

}

void SIRS::specMC_VD(int MC_cyc, double e_, double d_, double dI){ //MC with e and d
  MC_cycles = MC_cyc;
  d = d_;
  e = e_;
  d_I = dI;
  useVD = true;
}

/*
vec SIRS::derivativesSimple(vec yt){
  // S, I, R = y(0), y(1), y(2)
  dy(0) = c*(N-y(0)-y(1))
  dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  return dy;
}*/

vec SIRS::derivatives(vec yt){
  // S, I, R = y(0), y(1), y(2)
  dy(0) = c*(yt(2)) - a*yt(1)*yt(0)/N;
  dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  return dy;
}

vec SIRS::derivatives2(vec yt){
  //Update N since people die and get born
  // S, I, R = y(0), y(1), y(2)
  dy(0) = c*yt(2) - yt(0)*(a*yt(1)/N+d)+e*N;
  dy(1) = yt(1)*(a*yt(0)/N - b-d-d_I);
  dy(2) = b*yt(1) - (c+d)*yt(2);
  return dy;
}

vec SIRS::derivatives3(vec yt){
  // S, I, R = y(0), y(1), y(2)
  dy(0) = c*(yt(2)) - a*yt(1)*yt(0)/N - f;
  dy(1) = a*yt(0)*yt(1)/N - b*yt(1);
  dy(2) = b*yt(1) - c*yt(2) + f;
  return dy;
}

void SIRS::solveRK4(string filename){
  ofstream ofile;
  ofile.open(filename + "_RK4.csv");
  string prob_type = "std";
  if(useVD) prob_type = "VD";
  if(useV) prob_type = "Vac";
  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", RK4, " << prob_type << endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
  //for (double i = 0; i < t; i += dt){
  //double A = 1.5 ;double A0 = 4;double w = 0.08*(2*PI);
  for(int i = 0; i < num_pts; i++){
    //cout << i << endl;
    //a = A*cos(w*dt*i) + A0; //OPG D - SEASONAL VARIATION
    rk4(useVD);
    ofile << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
    //N = y(0) + y(1) + y(2);
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

  deadS = vec(num_pts, fill::zeros);
  deadI = vec(num_pts, fill::zeros);
  deadR = vec(num_pts, fill::zeros);
  deadI_dis = vec(num_pts, fill::zeros);
  S_born = vec(num_pts, fill::zeros);

  deadDis = vec(num_pts, fill::zeros);
  deadPop = vec(num_pts, fill::zeros);

  string prob_type;
  if(useVD) prob_type = "VD";
  else if(useV) prob_type = "Vac";
  else prob_type = "std";
  cout << "prob: " << prob_type << endl;

  ofile << t << ", " << dt << ", " << a << ", " << b << ", " << c << ", MC, " << prob_type  << endl;
  ofile << y(0) << ", " << y(1) << ", " <<  y(2);

  //Write death statistics to file
  if(useVD){
    ofile << deadDis(0) << ", " << deadPop(0) << ", " << S_born(0);
  }
  ofile << endl;
  cout << "dt: " << dt << endl;

  //double A = 1.5 ;double A0 = 4;double w = 0.08*(2*PI);
  double progress = 0;
  cout << num_pts << endl;
  for (int j = 0; j < MC_cycles; j ++){
    if(j+1 >= progress){
      cout << progress*100/MC_cycles << "% done." << endl; //Interval time: " << endl;
      progress += 0.1*MC_cycles;
    }
    int deadPopcount = 0;
    int deadDiscount = 0;
    int borncount = 0;
    for (int i = 1; i < num_pts; i++){
      //a = A*cos(w*dt*i) + A0; //OPG D - SEASONAL VARIATION
      //if(useSV) a = A*cos(w*dt*i) + A0; //OPG D - SEASONAL VARIATION
      MonteCarlo();
      deadPopcount += diedS + diedI + diedR;
      deadDiscount += diedI_disease;
      borncount += bornS;
      S_mc(i) += y(0);
      I_mc(i) += y(1);
      R_mc(i) += y(2);

      deadDis(i) += deadDiscount;
      deadPop(i) += deadPopcount;
      S_born(i) += borncount;
      //deadS(i) += deadS(i-1) + diedS;
      //deadI(i) += deadI(i-1) + diedI;
      //deadR(i) += deadR(i-1) + diedR;
      //deadI_dis(i) += deadI_dis(i-1) + diedI_disease;
      //S_born(i) += S_born(i-1) + bornS;
    }
    reset_states();
  }
  S_mc = S_mc/MC_cycles;
  I_mc = I_mc/MC_cycles;
  R_mc = R_mc/MC_cycles;

//dt = 1 dag. Simulere 1 år.
//På 1 år skal 7.8/1000 dø. Altså
//


  //deadS = deadS/MC_cycles;
  //deadI = deadI/MC_cycles;
  //deadR = deadR/MC_cycles;
  //deadI_dis = deadI_dis/MC_cycles;
  deadDis /= MC_cycles;
  deadPop /= MC_cycles;
  S_born = S_born/MC_cycles;

  cout << "Start write to file" << endl;
  for (int i = 1; i < num_pts; i++){
    ofile << S_mc(i) << ", " << I_mc(i) << ", " << R_mc(i);
    if(useVD){
      ofile << ", " << deadDis(i) << ", " << deadPop(i) << ", " << S_born(i);
    }
    ofile << endl;
  }
  cout << "done" << endl;
}


void SIRS::rk4(bool useVD){
  vec K1(3), K2(3), K3(3), K4(3);

  if(useVD == false && useV == false){
    //cout << "Using derivate1 function." << endl;
    K1 = dt*derivatives(y);
    //cout << "K1: " << K1 << endl;
    K2 = dt*derivatives(y+0.5*K1);
    K3 = dt*derivatives(y+0.5*K2);
    K4 = dt*derivatives(y+K3);
    y = y + (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0;
    y(2) = N - y(1) - y(0);
  }

  else if(useVD == true){
    //cout << "Using der2" << endl;
    K1 = dt*derivatives2(y);
    //cout << "K1: " << K1 << endl;
    K2 = dt*derivatives2(y+0.5*K1);
    K3 = dt*derivatives2(y+0.5*K2);
    K4 = dt*derivatives2(y+K3);
    y = y + (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0;
    N = y(0) + y(1) + y(2);
  }

  else if(useV == true){
    K1 = dt*derivatives3(y);
    K2 = dt*derivatives3(y+0.5*K1);
    K3 = dt*derivatives3(y+0.5*K2);
    K4 = dt*derivatives3(y+K3);
    y = y + (K1 + 2.0*K2 + 2.0*K3 + K4)/6.0;
  }
}

void SIRS::MonteCarlo(){
  // S, I, R = y(0), y(1), y(2)atom://teletype/portal/db276bfd-41b3-422c-bf85-e19906e7780d
  //Tansition probabilities and dt
  pR_S = (c*y(2))*dt;
  pS_I = (a*y(0)*y(1)/N)*dt;
  pI_R = (b*y(1))*dt;

  RS_count = SI_count = IR_count = 0;
  bornS = diedS = diedI = diedI_disease = diedR = 0;

  r = rand() % 100001;
  if (r/100000 < pR_S)
      RS_count ++;

  r = rand() % 100001;
  if (r/100000 < pS_I)
      SI_count ++;

  r = rand() % 100001;
  if (r/100000 < pI_R)
      IR_count ++;

  if(useVD){
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
    if (r/100000 < d_I*y(1)*dt) //Is birth rate given for dt = 1??
        diedI_disease ++;

    r = rand() % 100001;
    if (r/100000 < d*y(2)*dt) //Is birth rate given for dt = 1??
        diedR ++;

  }

  if(useV){
    r = rand() % 100001;
    if (r/100000 < f*y(0)*dt){
        y(2) ++;
        y(0) --;
    }
  }
  y(0) += RS_count - SI_count + bornS - diedS;
  y(1) += SI_count - IR_count - diedI - diedI_disease;
  y(2) += IR_count - RS_count - diedR;
  N = y(0) + y(1) + y(2); // Update tot pop after deaths and births
  //cout << "MC" << endl;
  //cout <<"S: " << y(0) << ", I: "<< y(1) << endl;
}




// inline void writeResults(fstream file);
//   file << y(0) << ", " << y(1) << ", " <<  y(2) << endl;
