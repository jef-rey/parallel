/* 
 * Iterative solutions to Linear Systems
 * Project 2 - CSCI 551
 * @author Jeff Anderson
 * @date 2/19/2020
 */

#define _XOPEN_SOURCE
/**
 * for some reason that i still haven't found an adequate answer to (that
 * i understand), i kept getting 'implicit declaration of function' errors for
 * srand48() and drand48(). looking into it, i found this #define _XOPEN_SOURCE
 * that seemed to resolve the situation. so i'm using it and its working and we
 * can all be happy about that.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h> // for srand

/*
 * @func: diag_dom
 * @desc: determines whether an n x n matrix is diagonally dominant or not
 * @return: a flag stating if the matrix is diagonally dominant or not
 */
int diag_dom(int n, double A[n][n]){
  // diagonal domiance means that x_ii would be greater than
  // all the other elements in that row

  int is_dominant = 1;// hoping for diagnoal dominance
  for(int i = 0; is_dominant && i < n; i++){
    // if we find a row that is not diagonally dominant, dump out immediately.
    double sum = 0;
    for(int j = 0; j < n; j++){
      if (i != j){
        sum += A[i][j];
      }
    }
    is_dominant = is_dominant && (sum < A[i][i]);
  }

  if(is_dominant){
    printf("the system is diagonally dominant\n");
  }else{
    printf("the system is not diagonally dominant\n");
  }

  return is_dominant;
}

double arae(int n,  double x_last[n], double x[n]);

/*
 * @func: jacobi
 * @desc: the jacobi method is a method to solve systems of linear equations. 
 *      it approximates a new estiamtion for each value x_1, x_2, ..., x_n
 *      and then uses those new estimations for new iterations until
 *      convergence is found
 *      my code was translated from the wikipedia page on Jacobi iterations,
 *      https://en.wikipedia.org/wiki/Jacobi_method
 */
void jacobi(int n, double A[n][n], double b[n], double x[n]){

  printf("running Jacobi iterations\n");
  int converge = 0; // flag that assumes no convergence
  double x_last[n]; // the previous values of x_1, x_2 .. x_n
  int e = 0; // declaring i here so i can use it as a counter
  for(; !converge && e < 500; e++){
    //check for convergence
    for(int i = 0; i < n; i++){
      x_last[i] = x[i];
      // mapping x_last to x = mapping current to previous (for first guess)
    } 
    // ACTUALLY CALCULATING JACOBI ITERATION
    for (int i = 0; i < n; i++){
      double sum = 0;
      for (int j = 0; j < n; j++){
        if (j != i){
          sum = sum + A[i][j] * x_last[j];
        }//if
      }//for
      x[i] = 1.0/A[i][i] * (b[i] - sum);

    }//for

    converge = converge || (arae(n, x_last, x) < 0.005);

  }
  if(converge){
    printf("system converged\n");
    printf("number of jacobi iterations = %d \n", e);

  }else{
    printf("system did not converge \n");
  }

  if (n < 11){
    for (int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        printf("%0.4lf ", A[i][j]);
      }
      printf("%0.4lf\n", b[i]);
    }

    for(int i = 0; i < n; i++){
      printf("%0.4lf ", x[i]);
    }
    printf("\n");
  }
}


/*
 * @func: gauss_seidel
 * @desc: Performs the Gauss Seidel method on a system of linear equations
 *    translated from the wikipedia article
 *    https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
 *    
 *    the Gauss Seidel method is extremely similar to the Jacobi method, except
 *    that it uses the new estimation of x_1 immediately in the new function of 
 *    x_2 instead of waiting for all iteration to be done 
 */
void gauss_seidel(int n, double A[n][n], double b[n], double x[n]){
  //iterate 500 times and check ARAE each time
  //if it ever hits 500, dump out and stop and report no convergence

  printf("running Gauss-Seidel iterations\n");
  int converge = 0; // flag that assumes no convergence
  double x_last[n]; // the previous values of x_1, x_2 .. x_n
  int e = 0; // declaring i here so i can use it as a counter
  for(; !converge && e < 500; e++){
    //check for convergence

    for (int i = 0; i < n; i++){
      x_last[i] = x[i];
      // mapping x_last to x = mapping current to previous (for first guess)
    }
    // ACTUALLY CALCULATING GAUSS-SEIDEL ITERATION
    for (int i = 0; i < n; i++){
      double sum = 0;
      for (int j = 0; j < n; j++){
        if (j != i){
          sum = sum + A[i][j] * x[j];
        }
      }
      x[i] = 1/A[i][i]* (b[i] - sum);
    }


    converge = converge || (arae(n, x_last, x) < 0.005);

  }
  if(converge){
    printf("system converged\n");
    printf("number of gauss_seidel iterations = %d \n", e);

  }else{
    printf("system did not converge \n");
  }

  if (n < 11){
    for (int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        printf("%0.4lf ", A[i][j]);
      }
      printf("%0.4lf\n", b[i]);
    }

    for(int i = 0; i < n; i++){
      printf("%0.4lf ", x[i]);
    }
    printf("\n");
  }
}

/*
 * @func: arae
 * @desc: finds the absolute relative approximate error
 * @returns the absolute relative approximate error
 */
double arae(int n, double x_last[n], double x[n]){
  int result = 0;
  for (int i = 0; i < n; i++){
    int candidate = ((x[i] - x_last[i]) / x[i] ); // * 100
    //candidate for ARAE = (curr - prev) / curr
    if (candidate < 0){
      candidate = candidate * -1; // abs 
    }
    if (candidate > result){
      result = candidate; // max -- from videos by mr. numerical guy
    }
  }
  return result;
}

/**
 * main driver of program
 */
int main(){

  int n; // size for arrays

  printf("read in the size of the linear system: ");
  scanf("%d", &n);

  srand48(time(0)); // seeding for drand48

  double A[n][n];
  double b[n];
  double x_jacobi[n];
  /** I had a single x array, but my Gauss Seidel was finishing too fast
   * because it already had the answer from running Jacobi. So I created
   * separate arrays to use to get more accurate results from solving a random
   * matrix of nxn size
   */
  double x_gs[n];

  //initialize vector x to all 1s for starting approximation
  for (int i = 0; i < n; i++){
    x_jacobi[i] = 1.0;
    x_gs[i] = 1.0;
  }

  if(n < 11){
    // read the augmented matrix for A | b from stdio, one row per line
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n+1; j++){

        double number;
        scanf("%lf", &number);
        if (j < n){
          A[i][j] = number;
        }
        b[i] = number;
      }
    }// A array values


  }else {
    /**fill A | b with random numbers using drand48 random num generator
     * scale result to the range of [-1.0e6, 1.0e6) if on the main diagonal
     * scale result to [-1.03e, 1.03e) if not on the diagonal.
     */
    for(int i = 0; i < n; i ++){
      for (int j = 0; j < n; j++){
        double number = drand48();
        number -= 0.5; 
        if (j == i){ // for the diagonal
          number = number * 2000000;
        }else{
          number = number * 2000;
        }
        A[i][j] = number;
      }
    } // filling A array with random values

    for(int i = 0; i < n; i++){
      double number = drand48();
      number -= 0.5; 
      number = number * 2000;
      b[i] = number;
    } // filling b array with random values

  }

  diag_dom(n, A);

  printf("\n");
  jacobi(n, A, b, x_jacobi);
  printf("\n");
  gauss_seidel(n, A, b, x_gs);

}

