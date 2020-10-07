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
    int N;
    double t_end;

    double totalKE;
    double totalPE;
    double pi = 2*acos(0.0);
    double G_scale = 4*pi*pi;

  public:

    //void addPlanet(Planet planet);
    void VerlocityVerlet();
    void EulerCromer();
    vec TotalAccelerationOnPlanet(int N_val, const Planet& planet, int index);
};

#endif
