#include "isingModel.hpp"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <functional>
#include <random>
#include <ostream>
using namespace std;

IsingModel::IsingModel(int n, double temp, int initMethod, long int numMC_cyc, double thread_seed){
  N = n;
  T0 = temp;
  M = 0;
  E = 0;
  numFlips = 0;
  sigma = 0;
  numMC_cycles = numMC_cyc;

  mt.seed(thread_seed);
  uniform_real_distribution<double> ddist(0,1);
  uniform_int_distribution<int> idist(0,N-1);
  //Set up spin_matrix NxN matrix, with spin up or down each element
  spin_matrix = new int*[N];
  for (int i = 0; i < N; i++)
      spin_matrix[i] = new int[N];

  expVals = new double[17]; //contains the possible deltaE vals
  for(int i=-8; i<=8; i+=4 ){
    expVals[i+8] = exp(-i/(T0));
  }
  //Initialze
  average = new double[5];
  for(int i = 0; i < 5; i++){
    average[i] = 0;
  }

  //Set up plus1 and min1 when using periodic boundary conditions
  plus1 = new int[N];
  min1 = new int[N];

  for(int i = 0; i < N; i++){
    plus1[i] = i+1;
    min1[i] = i-1;
  }
  plus1[N-1] = 0;
  min1[0] = N-1;

  //Fill initial matrix
  for (int i = 0; i < N; i ++){ // row
    for (int j = 0; j< N; j++ ){ //column
      if(initMethod == 0){
        spin_matrix[i][j] = 1;
      }
      if(initMethod == 1){
        spin_matrix[i][j] = -1;
      }

      if(initMethod == 2){
        double choice = ddist(mt);
        // 50-50 up or down
        if(choice < 0.5){
          spin_matrix[i][j] = -1;
        }
        else{
          spin_matrix[i][j] = 1;
        }
      }

      M += spin_matrix[i][j];
    }
  }
  findTotalEnergy(); // writes energy to E, class variable
  return;
}

void IsingModel::Metropolis(uniform_int_distribution<int> idist,uniform_real_distribution<double> ddist){

  int i_ = idist(mt);
  int j_ = idist(mt);
  double r = ddist(mt);
  double deltaE = calcE_ij(i_, j_);
  double e_exp = expVals[(int)deltaE+8];

  if(deltaE < 0){
    spin_matrix[i_][j_] *= -1;
    M = M + 2*spin_matrix[i_][j_];
    E = E + deltaE;
    numFlips += 1;
  }
  else if(r <= e_exp){
    spin_matrix[i_][j_] *= -1;
    M = M + 2*spin_matrix[i_][j_];
    E = E + deltaE;
    numFlips += 1;
  }

}

void IsingModel::solve(){
  //Choose random i and j and calculate the shift in E.
  //Calculate deltaE, if the random number r is <= exp(E/kT) then it happens
  //Confirm the flip and update spin matrix.
  double r;
  int N_sq = N*N;
  long int sampleCount = 0;
  ofstream ofile;
  double cutoff = 0.1;
  double loopCutoff = N_sq*cutoff*numMC_cycles;
  uniform_real_distribution<double> ddist(0,1);
  uniform_int_distribution<int> idist(0,N-1);

  for(long int i = 1; i <= loopCutoff; i++){
    Metropolis(idist, ddist);
  }
  for(long int i = loopCutoff; i <= (long int) N_sq*numMC_cycles ; i++){
    /*
    if(i > k*numMC_cycles*N_sq){
      cout << "Finish " << k*100 << " %" << endl;
      k += 0.1;
    }*/
    Metropolis(idist, ddist);
    sampleCount ++;
    average[0] += E; average[1] += E*E;
    average[2] += M; average[3] += M*M; average[4] += fabs(M);
  }

  for(int i = 0; i < 5; i++){
    average[i] /= (sampleCount); //divide value by the num of samples done in total
  }
  //Update important values when finished solving
  variance = (average[1] - average[0] * average[0])/N_sq;
  Cv = (average[1] - average[0] * average[0])/ (T0*T0) /N_sq;
  chi = (average[3] - average[4] * average[4]) / T0 / N_sq;
}

void IsingModel::writeFile(){
  ofstream Lfile;
  Lfile.open("Observables_" + to_string(N) + ".csv", ios_base::app);
  // T, <E>, <M>, Cv, chi
  Lfile << T0 << ", " << average[0] << ", " << average[4] << ", " << Cv << ", " << chi << endl;
}

void IsingModel::printValues(){
  //Print most important values to terminal
  cout << "Cv=" << Cv << ",  chi=" << chi << ",  variance=" << variance << endl;
  return;
}


double IsingModel::calcE_ij(int i, int j){
  // Calculate energy for group of 5 when potentially flipped.
  double E_now = 0;
  double E_after = 0;
  int s_ij = spin_matrix[i][j];
  int s_iplus = spin_matrix[plus1[i]][j];
  int s_jplus = spin_matrix[i][plus1[j]];
  int s_imin = spin_matrix[min1[i]][j];
  int s_jmin = spin_matrix[i][min1[j]];
  E_now -= s_ij*(s_iplus+s_jplus+s_imin+s_jmin);
  E_after = -E_now; // aka *=-1
  double deltaE = E_after-E_now; // -8,-4,0,4,8 possible values
  return deltaE;
}

void IsingModel::findTotalEnergy(){
  for (int i =0; i < N; i++){
    for (int j = 0; j < N; j++){
      double E_ij = spin_matrix[i][j]*spin_matrix[plus1[i]][j]+spin_matrix[i][j]*spin_matrix[i][plus1[j]];

      /* Print test for the 2x2 example
      cout << "Row +1: " <<  spin_matrix[plus1[i]][j] << endl;
      cout << "Col + 1: " << spin_matrix[i][plus1[j]] << endl;
      cout << "Energy for: i, j: " << i << ", " << j << "  . E: " << E_ij << "\n" << endl;
      */
      E -= E_ij; //energy for one point
    }
  }
}

IsingModel::~IsingModel(){
  //Frees up space when object is done.
  for(int i=0; i<N; i++){
    delete [] spin_matrix[i];
  }
  delete [] spin_matrix;
}

void IsingModel::printMatrix(){
  //Prints out full matrix to terminal
  cout << "Printed matrix" << endl;
  for (int i = 0; i < N; i ++){ // row
    for (int j = 0; j< N; j++ ){
      cout << spin_matrix[i][j] << " ";
    }
    cout << endl;
  }
}


void IsingModel::solve_write(){
  // Identical to solve, except writes to file for each sample to file
  double r;
  int N_sq = N*N;
  long int sampleCount = 0;
  ofstream ofile;
  double cutoff = 0.1;
  double loopCutoff = N_sq*cutoff*numMC_cycles;
  uniform_real_distribution<double> ddist(0,1);
  uniform_int_distribution<int> idist(0,N-1);

  ofile.open("e_hist.csv");
  long double k = 0.00; //Used if wish to see progress, used to get
  for(long int i = 1; i <= loopCutoff; i++){
    Metropolis(idist, ddist);
  }
  for(long int i = loopCutoff; i <= (long int) N_sq*numMC_cycles ; i++){
    Metropolis(idist, ddist);
    sampleCount ++;
    average[0] += E; average[1] += E*E;
    average[2] += M; average[3] += M*M; average[4] += fabs(M);
    //Write to file for each sample after cutoff
    ofile << average[0]/(sampleCount*N_sq) << ", " << average[4]/(sampleCount*N_sq) << ", " << numFlips << ", " << E/N_sq << endl;
  }

  for(int i = 0; i < 5; i++){
    average[i] /= (sampleCount); //divide value by the num of samples done in total
  }
  //Update important values when finished solving
  variance = (average[1] - average[0] * average[0])/N_sq;
  Cv = (average[1] - average[0] * average[0])/ (T0*T0) /N_sq;
  chi = (average[3] - average[4] * average[4]) / T0 / N_sq;
}
