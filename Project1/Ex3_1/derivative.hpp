#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP

#include <vector>
using namespace std;

class derivativeAtPoint {
  private:
    double m_a;
    vector<double> m_h, res2pt, res3pt;
    vector<float> m_h_f, res2pt_f, res3pt_f;

  public:
    void Initialize(vector<double> h,vector<float> h_f, double a);
    void FwEuler(double f(double x));
    void BackAndFw(double f(double x));
    void PrintResult();
};

#endif



/*
. Make thereafter a program which computes the first derivative using Eqs. (3.14) and (3.15)
as function of various step lengths h and let h → 0. Compare with the exact answer.
Your program should contain the following elements:
• A vector (array) which contains the step lengths. Use dynamic memory allocation.
• Vectors for the computed derivatives of Eqs. (3.14) and (3.15) for both single and double
precision.
• A function which computes the derivative and contains call by value and reference (for
C++ users only).
• Add a function which writes the results to file.
*/
