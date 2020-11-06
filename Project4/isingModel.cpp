#include "isingModel.hpp"
#include <iostream>
#include <math.h>
using namespace std;

void IsingModel::findTotalEnergy(){
  //int**& sp = spin_matrix;
  E = 0;
  cout << "Matrix: (2x2)" << endl;
  cout << spin_matrix[0][0] << " " << spin_matrix[0][1] <<  endl;
  cout << spin_matrix[1][0] << " " << spin_matrix[1][1] <<  endl;

  for (int i =0; i < N; i++){
    for (int j = 0; j < N; j++){
      // this is not right
      double E_ij = spin_matrix[i][j]+spin_matrix[plus1[i]][j]+spin_matrix[i][plus1[i]];
      cout << E_ij << endl;
      E -= E_ij; //energy for one point
    }
  }
  cout << "Energy in end: " << E << endl;
}

IsingModel::IsingModel(int n, double temp, int initMethod){
  N = n;
  T0 = temp;

  //Set up spin_matrix NxN matrix, with spin up or down each element
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
        spin_matrix[i][j] = rand() % 2;
      }

      M += spin_matrix[i][j];
    }
  }

  //Quick check if matrix set up right
  /*for (int i = 0; i < N; i ++){ // row
    for (int j = 0; j< N; j++ ){
      cout << i << ", " << j << "  Value: " << spin_matrix[i][j] << endl;
    }
  }*/

  //double energy = findTotalEnergy();
  return;
}

void IsingModel::solve(){

}

IsingModel::~IsingModel(){
  //Frees up space.
  for(int i=0; i<N; i++){
    delete [] spin_matrix[i];
  }
  delete [] spin_matrix;
}
