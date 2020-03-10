/**
 * program 3
 * parallel numerial integration with open mpi
 * Jeff Anderson
 * March 3, 2020
 *
 */

#include <stdlib.h>
#include <math.h>
// cosl, sinl, sqrtl, fabsl and powl for long double exp
#include <stdio.h>

/**
 * @func: f
 * @desc: f is the FUNCTION that we are hard coding to solve for this program,
 * tmin. it was provided in the writeup by Dr. Judy Challinger
 * @param: x - as per the equation, we are solving for x
 * @return: the value of the function at x
 */
long double f(long double x){
  
  return (-6*cosl(x/9) +sqrtl(powl((5*sinl(x/3) + 8*sinl(x/2)), 4)) + 10);

}

long double sum(long double a, long double b, int n){

  long double h;
  h = (b - a) / n; // taken directly from book
  long double area=(h/2)*(f(local_a) + f(local_b));// area of a single trapezoid
  long double local_a, local_b;

  long double sum = (f(a) + f(b)) / 2;
  // adapted from pacheo book
  for (int i = 1; i < n; i ++){
    local_a = a + i * h;
    sum += f(local_a);
  }
  sum *= h;

  return sum;
}

int main(){
  long double a = 0;
  long double b = 120;
  int n = 1; // init with 1 trapezoid at first RE: fig 3.4 from pacheo
  h = (b - a) / n; // taken directly from book
  long double area=(h/2)*(f(local_a) + f(local_b));// area of a single trapezoid

  printf("%LF\n", f(120));

  // use binary search for random guesses, since brute force is probably a bad
  // guess to start with (being that we were told that in the write up. )
  
  long double stop_criteria = 0.5E-14;

  //TODO: make loop for stopping criteria

  


}
