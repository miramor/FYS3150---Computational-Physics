#include "derivative.hpp"
#include <iostream>
#include <cmath>

using namespace std;

void derivativeAtPoint::Initialize(vector<double> h, vector<float> h_f, double a){
  //Input for h values: {a,b,c,d,e}
  m_h = h;
  m_h_f = h_f;
  m_a = a;
}

void derivativeAtPoint::FwEuler(double f(double x)){
  double a = m_a;
  vector<double> h = m_h;
  vector<float> h_f = m_h_f;
  //h, h_f = m_a, m_h, m_h_f;
  for(int i = 0; i < h.size(); i++){
    double h_val = h.at(i);
    double res_i = (f(a+h_val)-f(a))/h_val;
    res2pt.push_back(res_i);

    float h_fval = h_f.at(i);
    double res_i_f = (f(a+h_fval)-f(a))/h_fval;
    res2pt_f.push_back(res_i_f);
  }
}

void derivativeAtPoint::BackAndFw(double f(double x)){
  //a, h, h_f = m_a, m_h, m_h_f;
  double a = m_a;
  vector<double> h = m_h;
  vector<float> h_f = m_h_f;

  for(int i = 0; i < m_h.size(); i++){
    double h_val = m_h.at(i);
    double res_i = (f(a+h_val)-f(a))/h_val;
    res3pt.push_back(res_i);

    float h_fval = m_h_f.at(i);
    double res_i_f = (f(m_a+h_fval)-f(m_a))/h_fval;
    res3pt_f.push_back(res_i_f);
  }
}

void derivativeAtPoint::PrintResult(){
  cout << "hi" << endl;

  for (auto i = res2pt.begin(); i != res2pt.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  for (auto i = res3pt.begin(); i != res3pt.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  for (auto i = res2pt_f.begin(); i != res2pt_f.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  for (auto i = res3pt_f.begin(); i != res3pt_f.end(); ++i)
    cout << *i << ' ';
  cout << endl;

  /*
  cout << "Choose method: (1) FwEuler, (2)BackAndFw, (3) Both";
  int cin >> choice;
  if (choice == 1){
    cout << "Result for FwEul" << res2pt << endl;
    cout << "Float: " << res2pt_f << endl;
  }
  else if (choice == 2){
    cout << "Result for BackAndFw" << res3pt << endl;
    cout << "Float: " << res3pt_f << endl;
  }
  else if (choice == 3){
    cout << "Result for FwEul" << res2pt << endl;
    cout << "Result for BackAndFw" << res3pt << endl;
  }
  else {
    cout << "Bad input, please put 1 or 2 instead of: " << choice << endl;
  }*/
}
