using namespace std;
# include <iostream>
// begin main functio

void func(int x, int *y);

int main(int argc, char argv[])
  {
    int a;
    int *b;
    a = 10;
    b = new int[10];

    for( int i = 0; i < 10; i++){
      b[i] = i;
  }
  cout << b << endl;
  cout << *b << endl;
  cout << b[6] << endl;

  func(a,b);
  return 0;
} // end of main function
  // definition of the function func

void func(int x, int *y){
  x += 7;
  *y += 10;
  y[6] += 10;

  cout << y << endl;
  cout << *y << endl;
  cout << y[6] << endl;
  return;
} // end function func
