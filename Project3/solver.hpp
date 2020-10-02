#ifndef PLANET.HPP
#define PLANET_HPP

#include <vector>
#include <string.h>
#include <armadillo>

class Solver{
  private:

    friend class planet;

    int totalPlanets;
    vec<planet> planets;

    double totalKE;
    double totalPE;
  public:

    addPlanet(Planet planet);
}

#endif
