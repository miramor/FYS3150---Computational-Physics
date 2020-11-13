#include "isingModel.hpp"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
using namespace std;

void IsingModel::findTotalEnergy(){
  for (int i =0; i < N; i++){
    for (int j = 0; j < N; j++){
      // this is not right

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

IsingModel::IsingModel(int n, double temp, int initMethod){
  N = n;
  T0 = temp;
  M = 0;
  E = 0;
  numFlips = 0;

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
  srand (time(NULL)); // Set seed for random gen numbers

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
        double choice = drand48();

        //cout << choice << endl;
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

void IsingModel::Metropolis(){
  double a =  drand48();
  double b = drand48();
  int i_ = (int)(a*N);
  int j_ = (int)(b*N);

  //int i_ = rand() % N; //error when defined i_ outside loop
  //int j_ = rand() % N;
  double r = drand48();
  double deltaE = calcE_ij(i_, j_);
  double e_exp = expVals[(int)deltaE+8];

  if(deltaE < 0){
    //cout << "Flip success!" << endl;
    spin_matrix[i_][j_] *= -1;
    M = M + 2*spin_matrix[i_][j_];
    E = E + deltaE;
    numFlips += 1;
  }
  //cout << "r: " << r << "  , exp: " << e_exp << " E: " << deltaE << endl;
  else if(r <= e_exp){
    //cout << "Flip success!" << endl;
    spin_matrix[i_][j_] *= -1;
    M = M + 2*spin_matrix[i_][j_];
    E = E + deltaE;
    numFlips += 1;
  }
}

void IsingModel::solve(){
  srand (time(NULL)); // Set seed for random gen numbers
  //Choose random i and j and calculate the shift in E.
  // Calculate deltaE, if the random number r is <= exp(E/kT) then it happens
  // Confirm the flip and update spin matrix.
  double r;
  int N_sq = N*N;
  int numMC_cycles = 1000000;  // num of monte carco cycles
  int sampleCount = 0;
  //N_sq = 2;

  ofstream ofile;
  ofile.open("e_hist.csv");
  double cutoff = 0.1;
  double loopCutoff = N_sq*cutoff*numMC_cycles;
  ofile << cutoff << ", " << numMC_cycles << ", " << T0 << ", " << N << endl;

  for(int i = 1; i <= loopCutoff; i++){
    Metropolis();
  }

  for(int i = loopCutoff; i <= (long int) N_sq*numMC_cycles ; i++){
    Metropolis();
    sampleCount ++;
    average[0] += E; average[1] += E*E;
    average[2] += M; average[3] += M*M; average[4] += fabs(M);
    //ofile << average[0]/sampleCount << ", " << average[4]/sampleCount << ", " << numFlips << ", " << E/N_sq << endl;

  }
  for(int i = 0; i < 5; i++){
    average[i] /= (sampleCount);
  }

  Cv = (average[1] - average[0] * average[0]) / T0 / N_sq;
  chi = (average[3] - average[4] * average[4]) / (T0*T0) / N_sq;

  for(int i = 0; i < 5; i++){
    average[i] /= N_sq;
  }

  //cout << "<E> " << average[0]<< "  <M> " << average[4] << endl;
  //cout << "CV: " << Cv << "  chi: " << chi << endl;
}

void IsingModel::writeFile(){
  ofstream Lfile;
  Lfile.open("Observables_" + to_string(N) + ".csv", ios_base::app);
  // T, <E>, <M>, Cv, chi
  Lfile << T0 << ", " << average[0] << ", " << average[4] << ", " << Cv << ", " << chi << endl;
}



double IsingModel::calcE_ij(int i, int j){
  // check E when flipped, does not actually modify the spin_matrix
  double E_now = 0;
  double E_after = 0;
  int s_ij = spin_matrix[i][j];
  int s_iplus = spin_matrix[plus1[i]][j];
  int s_jplus = spin_matrix[i][plus1[j]];
  int s_imin = spin_matrix[min1[i]][j];
  int s_jmin = spin_matrix[i][min1[j]];

  E_now -= s_ij*(s_iplus+s_jplus+s_imin+s_jmin);
  E_after = -E_now; //tilsvarer å gange s_ij med -1

  double deltaE = E_after-E_now; // -8,-4,0,4,8 possible values
  return deltaE;
}

IsingModel::~IsingModel(){
  //Frees up space.
  for(int i=0; i<N; i++){
    delete [] spin_matrix[i];
  }
  delete [] spin_matrix;
  //cout << "Freed up space" << endl;
}

double* IsingModel::getAverage(){
  return average;
}

void IsingModel::printMatrix(){
  cout << "Printed matrix" << endl;
  for (int i = 0; i < N; i ++){ // row
    for (int j = 0; j< N; j++ ){
      cout << spin_matrix[i][j] << " ";
    }
    cout << endl;
  }
}
