#ifndef ISINGMODEL_HPP
#define ISINGMODEL_HPP

#include <vector>
#include <string.h>
#include <random>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <functional>
#include <ostream>
using namespace std;


class IsingModel{
  private:
    int** spin_matrix;//spin_matrix[i][j]
    int numFlips;

    double* E_array; //history of all energies
    double* deltaE_vals;
    double* average;
    double* expVals;

    int* plus1; //Used for indexing neighbours (i+1) and (j+1)
    int* min1;

    double E; // current energy in system
    double M; // magnetization
    double Cv; //heat capacity
    double chi; //Susceptibility
    double sigma; //standard deviation for energy
    double variance; //variance (sigma**2) for energies

    double T0; //temperatur, also equals T*k_b since k = 1
    int N; // Size of matrix (L)
    long int numMC_cycles; //Num cycles to execute

    mt19937 mt; // Genrator for random number generator

  public:
    IsingModel(int n, double temp, int initMethod, long int numMC_cyc, double thread_seed); // initMethod: (0)up, (1)down or (2)random
    void findNeighbour();
    void swapSpinOnce();
    void findTotalEnergy();
    void writeFile();
    void solve();
    void printMatrix();
    void Metropolis(uniform_int_distribution<int> idist, uniform_real_distribution<double> ddist);
    double* getAverage(); //contains <E>, <M>, <E**2> used for plotting
    double calcE_ij(int i, int j);
    void printValues();
    ~IsingModel(); //deletes array when object is "dead"
};

#endif
