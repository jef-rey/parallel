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

long double approx_area(long double a, long double b, int n){

  long double h = (b - a) / n; // taken directly from book
  // n is the number of trapezoids we're using at the particular iteration
  //long double area=(h/2)*(f(local_a) + f(local_b));// area of a single trapezoid
  long double local_a;

  long double sum = (f(a) + f(b)) / 2;
  // adapted from pacheo book
  for (int i = 1; i < n; i ++){
    local_a = a + i * h;
    //local_a is really like local_b...
    //in that local_b the next point to iterate on.
    sum += f(local_a);
  }
  sum *= h;

  return sum;
}

long double arte(long double approx_value){

  long double true_value = 6745.9948824470072783; // 20 digits
  long double true_error = true_value - approx_value;
  long double result = fabsl(true_error / true_value);

  return result;

}

void report(long double a, long double d, int t){

  long double approx = approx_area(a, b, t);
  long double abs_rel_tru_error = arte(approx); 
  printf("the integration result = %.13Le\n", approx);
  //USING Le INSTEAD OF e AS e ONLY RESULTS IN WARNINGS WHEN USING LONG DOUBLES
  printf("the absolute relative true error is = %.19Le\n", abs_rel_tru_error);
}

int find_UB(long double a, long double b, int t){
  long double stop_criteria = 0.5E-14;
  while(arte(approx_area(a, b , t) > stop_criteria)){
    // if this returns true, then we have _an_ upperbound.
    t *= 2;
    // doubling t to get that upper bound. 
  }

  return t;
}

void search(long double a, long double b, int t){

  int upper_bound = find_UB(a, b, t);
  int lower_bound = 1;

  if (t < 1){
    fprintf(stderr, "t cannot be < 1 ");
  }

}

int binary_search(long double a, long double b, int lower_bound, int upper_bound ){

  long double stop_criteria = 0.5E-14;
  int mid = lower_bound + (upper_bound - lower_bound ) / 2;

  if (upper_bound >= lower_bound){
    //do something
    //if the upper bound
    if (mid == x) 
      return mid; 

    //If element is smaller than mid, then 
    // it can only be present in left subarray 
    if (arte(approx_area(a,b,mid)) <= stop_criteria){ 
      // then mid is definitely not a valid answer
      return binary_search(a, b, lower_bound, mid); 
    }
    // else the element can only be present in the R subarrahy
      return binary_search(a, b, mid + 1, upper_bound);

  }



}

int main(){
  //h = (b - a) / n; // taken directly from book
  //long double area=(h/2)*(f(local_a) + f(local_b));// area of a single trapezoid

  printf("%LF\n", f(120));

  // use binary search for random guesses, since brute force is probably a bad
  // guess to start with (being that we were told that in the write up. )

  //long double stop_criteria = 0.5E-14;

  //TODO: make loop for stopping criteria
  //comment for testing
  //TODO: CLARITY ON WHY WE NEED TO READ IN A AND B FROM STDIN IF WE'RE
  //HARDCODING THIS FOR A PARTICULAR FORMULA.
  //anything other than a = 0, and b = 120 will result in a different
  //true_value
  //ALSO -- for [R]eporting, are we supposed to be printing the _absolute_
  //relative true error, or the relative true error, or the true error

  /*
   * absolute relative true error is:
   * abs( relative true error ) = true error / true value
   *
   */

  char input;
  long double a = 0;
  long double b = 120;
  int t = 1; // init with 1 trapezoid at first RE: fig 3.4 from pacheo
  printf("select [R]eport or [S]earch: ");
  scanf("%c", &input);
  if (input == 'R'){
    //reporting
    printf("type in values for integration: lower bound, upper bound and number of [T]raps: ");
    scanf("%LF %LF %d", &a, &b, &t);
    report(a, b, t);

  }else if (input == 'S'){
    //stuff for S
  }else{
    printf("I don't know what you mean.");
    exit (1);
  }





}


/** binary search notes
 * given a starting t value, if the approx using t has an ARTE <= stopping
 * criteria, then it's likely that t is probably too big. 
 * so then you need to look at the midpoint below t. mid = 1 - (t - 1) / 2
 * keep searching recursively like that.
 * once you go a step down, now you have a new inclusive upper bound.
 *
 * what if the original t is too small?
 * then you have to keep doubling everything, until you find an ARTE that is <=
 * stopping criteria. 
 * then you have a new upper bound, and a lower bound (the last thing you
 * tried). 
 * new bounds will be > old lower bound, might <= upper bound. 
 * at which point, it's the same search criteria as the case above. 
 *
 * ub = # of traps
 * 
 * is arte(approx_are(a, b, t)) <= stop_criteria
 */
