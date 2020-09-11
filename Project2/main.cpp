int main(int argc, char const *argv[]) {

  //Run a loop to rotate until we reach max iterations or reach limit
  while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
    max:offdiag = maxoffdiag ( A, &k, &l, n );
    rotate ( A, R, k, l, n );
    iterations++;
  return 0;
}
