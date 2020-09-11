#ifndef JACOBIALGO_HPP
#define JACOBIALGO_HPP

#include <vector>
#include <string.h>
#include <tuple>

class JacobiEigenSolve {
  private:
    Mat<double> A; // Symmetric matrix to diagonalize,  Dim: n x n
    Mat<double> V; //Matrix to contain eigenvectors
    int n; // Size of the matrix
    double eps = 1.0e-8;
    double max_iterations = (double) n * (double) n * (double) n;
    int iterations = 0;
    double a, b, c; // a - lower diagonal, b - middle diagonal, c - upper diagonal
    //int k, l; // Changed each time
  public:
    void Initialize(double a_val, double b_val, double c_val, int max_ite, int n_val); // make the symmetric matrix and empty V matrix.
    tuple<double, int, int> FindMaxEle(Mat<double> A); // Returns a tuple of value
    //of element and its index. (val, row, column) => f.eks (2.2, 1, 3)
    void Rotate(Mat<double> A, int row, int col); // finds the value of cos(theta) and sin(theta)
    void Solve(); //Simply runs the rotation until we reached the max ite or reached the eps
    Mat<double> PrintA(); // Prints A, used for checks

    //enkel matrise 3x3, analyisk rotasjonsmatrisen. Gjøre alle steg for hånd
    // sjekke om egenverdiene funker.
};

#endif
