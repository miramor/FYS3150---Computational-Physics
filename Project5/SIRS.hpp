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
    double d;
    double e;
    double d_I;
    int t;
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

    void rk4(vec (*)(vec));
    vec derivatives1(vec yt);
    vec derivatives2(vec yt);
    void MonteCarlo();
    void reset_states();


  public:
    SIRS(double S, double I, double a_, double b_, double c_, double t_, double dt_); //RK4
    SIRS(double S, double I, double a_, double b_, double c_, double t_, double dt_,
       double e_, double d_, double dI); //RK4 with e and d
    SIRS(double S, double I, double a_, double b_, double c_, double t_, int MC_cyc); //MC
    SIRS(double S, double I, double a_, double b_, double c_, double t_, int MC_cyc,
       double e_, double d_, double dI); //MC with e, d and d_I

    void solveMC(string filename);
    void solveRK4(string filename);
};




#endif
