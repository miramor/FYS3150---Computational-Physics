#ifndef ISINGMODEL_HPP
#define ISINGMODEL_HPP

#include <vector>
#include <string.h>
using namespace std;


class IsingModel{
  private:
    int** spin_matrix;//spin_matrix[i][j]
    int monteCycles; // nr of cycles the system has been through, N**2

    double* E_array; //history of all energies
    double* deltaE_vals;

    int* plus1;
    int* min1;

    double E; // current energy in system
    double E_exp;
    double E_sq_exp;
    double M; // magnetization
    double Cv; //heat capacity

    double T0; //temperatur, also equals T*k_b since k = 1
    double T_end;
    int N;

  public:
    IsingModel(int n, double temp, int initMethod); // initMethod: (0)up, (1)down or (2)random
    void findNeighbour();
    void swapSpinOnce();
    void findTotalEnergy();
    void solve();
    ~IsingModel(); //deletes array when object is "dead"
};

#endif
