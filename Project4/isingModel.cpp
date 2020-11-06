#include "isingModel.hpp"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

void IsingModel::findTotalEnergy(){
  //int**& sp = spin_matrix;
  E = 0;
  /*
  cout << "Matrix: (2x2)" << endl;
  cout << spin_matrix[0][0] << " " << spin_matrix[0][1] <<  endl;
  cout << spin_matrix[1][0] << " " << spin_matrix[1][1] <<  endl;
  */

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
  //cout << "Energy in end: " << E << endl;
}

IsingModel::IsingModel(int n, double temp, int initMethod){
  N = n;
  T0 = temp; //er du her?

  //Set up spin_matrix NxN matrix, with spin up or down each element
  spin_matrix = new int*[N];
  for (int i = 0; i < N; i++)
      spin_matrix[i] = new int[N];

  deltaE_vals = new double[5]; //contains the possible deltaE vals
  double E_vals[] = {-8,-4,0,4,8};

  for(int i=0; i<5; i++){
    deltaE_vals[i] = exp(E_vals[i]/(T0));
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
        int choice = rand() % 2;
        //cout << choice << endl;
        if(choice == 0){
          choice = -1;
        }
        spin_matrix[i][j] = choice;
      }

      M += spin_matrix[i][j];
    }
  }

  //Print out matrix for visual tests:
  findTotalEnergy(); // writes energy to E, class variable
  cout << "Matrix: " << "E: " << E << " , M: " << M << endl;
  for(int i = 0; i<N ; i++){
    for(int j = 0; j<N ; j++){
      cout << spin_matrix[i][j] << "  ";
    }
    cout << endl;
  }

  //Quick check if matrix set up right
  /*for (int i = 0; i < N; i ++){ // row
    for (int j = 0; j< N; j++ ){
      cout << i << ", " << j << "  Value: " << spin_matrix[i][j] << endl;
    }
  }*/

  return;
}

void IsingModel::solve(){
  srand (time(NULL)); // Set seed for random gen numbers
  //Choose random i and j and calculate the shift in E.
  // Calculate deltaE, if the random number r is <= exp(E/kT) then it happens
  // Confirm the flip and update spin matrix.
  int rint = rand() % 101;


  int i_ = rand() % N;
  int j_ = rand() % N;
  cout << "Random i and j:" << endl;
  cout << i_ << " " << j_ << endl;

  double deltaE = calcE_ij(i_, j_);

  double r = (double) rint/100;// between 0 and 1, 100 possible pts

  double e_exp = exp(deltaE/T0);
  if(r < e_exp){
    spin_matrix[i][j] *= -1;
  }
}

double IsingModel::calcE_ij(int i, int j){
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
  cout << "E0: " << E_now << endl;
  cout << "E_after: " << E_after << endl;
  cout << "deltaE: " << deltaE << endl;

  // nå enten hent den boltzmann verdien fra map eller array

  //temp sol




  return deltaE;
}

IsingModel::~IsingModel(){
  //Frees up space.
  for(int i=0; i<N; i++){
    delete [] spin_matrix[i];
  }
  delete [] spin_matrix;
}
