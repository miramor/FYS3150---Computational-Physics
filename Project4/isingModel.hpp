#ifndef ISINGMODEL_HPP
#define ISINGMODEL_HPP

#include <vector>
#include <string.h>
using namespace std;


class IsingModel{
  private:
    int** spin_matrix;//spin_matrix[i][j]
    int numFlips;

    double* E_array; //history of all energies
    double* deltaE_vals;
    double* average;
    double* expVals;

    int* plus1;
    int* min1;

    double E; // current energy in system
    double E_exp;
    double E_sq_exp;
    double M; // magnetization
    double Cv; //heat capacity
    double chi; //Susceptibility

    double T0; //temperatur, also equals T*k_b since k = 1
    double T_end;
    int N;

  public:
    IsingModel(int n, double temp, int initMethod); // initMethod: (0)up, (1)down or (2)random
    void findNeighbour();
    void swapSpinOnce();
    void findTotalEnergy();
    void writeFile();
    void solve();
    void printMatrix();
    void Metropolis();
    double* getAverage(); //contains <E>, <M>, <E**2> used for plotting
    double calcE_ij(int i, int j);
    ~IsingModel(); //deletes array when object is "dead"
};

#endif
