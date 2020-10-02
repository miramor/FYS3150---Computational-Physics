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
  public:

    //void addPlanet(Planet planet);
    void verletSolve();
}

#endif
