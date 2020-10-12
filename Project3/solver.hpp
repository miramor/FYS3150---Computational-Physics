#ifndef SOLVER_HPP
#define SOLVER_HPP


#include <string.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <vector>
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
    string sysName;

  public:
    Solver(vector<Planet> sysPlanets, int N_val, double t_n_val, string sys);
    //void addPlanet(Planet planet);
    void VelocityVerlet();
    void EulerCromer();
    void Euler();
    vec TotalAccelerationOnPlanet(Planet& planet, int index);
    void WriteResults();

    double getTotalEnergy();
};

#endif
