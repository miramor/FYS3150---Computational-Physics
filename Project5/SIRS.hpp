#ifndef SIRS_HPP
#define SIRS_HPP
#include <string>
#include <armadillo>
using namespace arma;
using namespace std;

class SIRS{
  private:
    vec y;
    double a;
    double b;
    double c;
    double t;
    double dt;
    int num_pts;
    double N;

    //MC variables
    int MC_cycles;
    vec S_mc;
    vec I_mc;
    vec R_mc;
    double pS_I, pI_R, pR_S;
    double r;

    int S0, I0, R0; //used to reset initial state when repeating MC

    //RK4 variables
    vec dy;
    bool useDer1;

    //VD variables
    double d;
    double e;
    double d_I;

    void rk4(bool useDer1);
    vec derivatives(vec yt);
    vec derivatives2(vec yt);
    void MonteCarlo(bool useVD);
    void reset_states();


  public:
    SIRS(double S_, double I_, double a_, double b_, double c_, double t_); //General

    void specRK4(double dt_); //RK4 std
    void specMC(int MC_cyc); // MC std

    void specRK4_VD(double dt_, double e_, double d_, double dI); //RK4 vital dynamics
    void specMC_VD(int MC_cyc, double e_, double d_, double dI); //MC vital dynamics

    void solveMC(string filename);
    void solveRK4(string filename);
};




#endif
