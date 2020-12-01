#ifndef SIRS_HPP
#define SIRS_HPP
#include <string>
#include <armadillo>
using namespace arma;
using namespace std;

class SIRS{
  private:
    //vector <double> y, dy;
    vec y;
    vec dy;
    double a;
    double b;
    double c;
    int t;
    double dt;
    double dt2;
    double N;


    void rk4();
    //vector <double> derivatives(vector <double> yt);
    vec derivatives(vec yt);


  public:
    SIRS(double S, double I, double a_, double b_, double c_, double t_, double dt_);
    void solve(string filename);
};




#endif
