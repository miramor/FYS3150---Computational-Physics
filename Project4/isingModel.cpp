int * make_spinMatrix(int N){
  int N_2 = N+2;
  int matrix[N_2*N_2];
  //int matrix[N*N];

  /*
  for (int i = 0; i < N*N; i ++){
    if (rand() % 2 == 1)
      matrix[i] = 1;
    else matrix[i] = -1;
  }
  */

  for (int i = 1; i < N_2-2; i ++){ // row
    for (int j = 1; j; N_2-2, j++ ){ //colum
      if (rand() % 2 == 1) matrix[i][j] = 1;
      else matrix[i][j] = -1;
      M += matrix[i][j];
    }
  }


  for (int k = 1; k < N_2-2; k++){
    matrix[0][k] = matrix[N_2-2][k];
    matrix[N_2-1][k] = matrix[1][k];
    matrix[k][N_2-1] = matrix[k][1];
    matrix[k][0] = matrix[k][N_2-2];
    }

    for (int i =1; i < N_2-2, i++){
      for (int j = 1, j< N_2, j++){
        E -= matrix[i][j]*( matrix[i][j+1] + matrix[i+1][j]) + matrix[i+1][j+1]*(matrix[i][j+1] + matrix[i+1][j]);
      }
    }
    E *= 2*J;
  cout << matrix[1] << matrix[2]<< endl;
  return matrix;
}

#include "IsingModel.hpp"
using namespace std;


IsingModel::IsingModel(int n, double temp, int initMethod){
  N = n;
  T0 = temp;

  spin_matrix = 
}

IsingModel::FindNeighbour(int i, int j){

}
