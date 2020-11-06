#ifndef PLANET_HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>
using namespace arma;
using namespace std;


class IsingModel{
  private:
    int** spin_matrix;//spin_matrix[i][j]
    int monteCycles;

    double[] E_array;
    double E;
    double E_exp;
    double E_sq_exp;
    double M;

    double T0;
    double T_end;
    int N;

  public:
    IsingModel(int n, double temp, int initMethod); // initMethod: (0)up, (1)down or (2)random
}

#endif
