#ifndef SOLVER_HPP
#define SOLVER_HPP


#include <string.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <vector>
#include "planet.hpp"

using std::vector;

class Solver{
  private:

    friend class Planet;

    int totalPlanets;
    vector<Planet> planets;

    double totalKE;
    double totalPE;
    double pi = 2*acos(0.0);
    double G_scale = 4*pi*pi;
  public:

    //void addPlanet(Planet planet);
    void verletSolve();
    void EulerCromer(int N, double t_end, ve<Planet> planets);
    double TotalForceOnPlanet(const Planet& planet);
};

#endif
