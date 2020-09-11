#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "jacobiAlgo.hpp"
#include "time.h"

using namespace std;
using namespace arma;


void JacobiEigenSolve::Initialize(double a_val, double b_val, double c_val, int max_ite, int n_val){

}

// Returns a tuple of value of element and its index. (val, row, column) => f.eks (2.2, 1, 3)
tuple<double, int, int> FindMaxEle(Mat<double> A){


}

// finds the value of cos(theta) and sin(theta) and fills up V each time
void Rotate(Mat<double> A, int row, int col){

}

 //Runs the rotation until we reached the max ite or reached the eps
void Solve(){

}

// Prints A, used for checks
Mat<double> PrintA(){
  cout << "Matrix A: " << A << endl;
}
