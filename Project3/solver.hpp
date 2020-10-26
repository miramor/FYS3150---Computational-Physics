#ifndef SOLVER_HPP
#define SOLVER_HPP


#include <string.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "planet.hpp"

using std::vector;
using namespace std;

class Solver{
  private:

    friend class Planet;

    int totalPlanets;
    vector<Planet> planets;
    int N;
    int twoN;
    double t_n;
    double totalKE;
    double totalPE;
    double pi = 2*acos(0.0);
    double G_scale = 4*pi*pi;
    double G = 6.67408e-11; //Gravitational constant
    string sysName;
    string method;

  public:
    Solver(vector<Planet> sysPlanets, int N_val, double t_n_val, string sys);
    //void addPlanet(Planet planet);
    void VelocityVerlet();
    void VelocityVerlet2();
    void EulerCromer();
    void Euler();
    vec TotalAccelerationOnPlanet(Planet& planet, int index);
    void WriteResults();
    void WritePeriResults();
    void testTotE();
    void testAngMom();
    double calcPE(int k, int j);
    double calcKE(int k, int j);
    double calcL(int index);

    void VertleNoStorage();
    vec TotalAccelerationOnPlanet_opt(Planet& planet, bool useCurr);

};

#endif
