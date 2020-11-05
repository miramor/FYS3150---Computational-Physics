#include <cstdio>
#include <omp.h>
#include <iostream>
#define NUM_THREADS 2
using namespace std;

int main(int argc, char const *argv[]) {
  static long num_steps = 10; double step;
  double pi = 0.0;
  step = 1.0/(double) num_steps;
  omp_set_num_threads(NUM_THREADS);
  int i, id, nthrds;
  double x, sum;

  #pragma omp parallel
  {
    id = omp_get_thread_num();
    cout << id << endl;
    nthrds = omp_get_thread_num();
    //if (id==0) nthreads = nthrds;
    for (i=id; i< num_steps; i=i+nthrds){
      sum=0.0;
      x = (1+0.5)*step;
      sum+= 4.0/(1.0+x*x);
    }
  }
    #pragma omp atomic
      pi += sum;
    cout << sum << endl;
  return 0;
}
//-----------------------------------
#pragma omp parallel for reduction (+:ave) //video 10
//where are barriers?
//work sharing construct
//Lock Routines
