#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>
#include <string.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "planet.hpp"
class Solver{
  private:

    friend class Planet;

    int totalPlanets;
    std::vector<Planet> planets;

    double totalKE;
    double totalPE;
  public:

    //void addPlanet(Planet planet);
    void verletSolve();
};

#endif
