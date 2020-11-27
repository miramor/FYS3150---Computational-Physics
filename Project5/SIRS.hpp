#ifndef SIRS_HPP
#define SIRS_HPP


using namespace std;

class SIRS{
  private:
    double S;
    double I;
    double R;

    double a;
    double b;
    double c;

  public:
  SIRS(double S, double I, double a, double b, double c);
  void rk4();
  void solve();
};




#endif
