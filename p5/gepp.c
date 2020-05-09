/* p5
* Jeff Anderson
* CSCI 551 - Gaussian Elimination with Partial Pivoting
* May 2020
*
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

/**
 * \brief L2 norm computes the square root of the sum of the squares of the
 * residual vector Ax-b
 * \params n, the size of the matrix
 *         A, the LHS of the matrix, filled with random variables
 *         b, the RHS of the matrix, here as a vector of answers
 *         x, the unknowns
 * \return the sum of the squares of the residual vector Ax-b
 */
double L2norm(int n, double* A, double* b, double* x){

  double Ax[n];

  for(int i = 0; i < n; i++){
    double sum = 0.0;
    for (int j=0; j < n; j++){
      sum += A[i*n+j] * x[j];

    }
    Ax[i] = sum;
  }

  double squares_sum = 0.0;
  for (int i = 0; i < n; i++){
    double difference= Ax[i] - b[i];
    squares_sum += difference * difference; 

  }
  return sqrt(squares_sum);
}

/**
 * \brief function to help with partial pivoting in gaussian elimination
 * \decription we use pivoting to swap the current row with the row with the
 * largest absolute value in the column of interest. this technique avoids
 * division by zero errors and reduces round off errors
 * \params n, the size of the matrix
 *         A, the LHS of the matrix, filled with random variables
 *         b, the RHS of the matrix, here as a vector of answers
 *         diag, the [i][j] position that we are looking at
 * \return void function
 */
void pivot(int n, double* A, double* b, int diag ){
  int i;

  double greatest_val = 0;
  int row_with_greatest_elem = diag;

  for (/*int*/ i = diag; i < n; i++){
    double abs_elem = fabs(A[i*n+diag]);
    if (abs_elem > greatest_val){
      greatest_val = abs_elem;
      row_with_greatest_elem = i;
    }
  }

  //swap rows in A column
  for (/*int*/ i = diag; i < n; i++){
    double tmp = A[diag*n+i];
    A[diag*n+i] = A[row_with_greatest_elem*n+i];
    A[row_with_greatest_elem*n+i] = tmp;
  }

  double temp;
  temp = b[diag] ;
  b[diag] = b[row_with_greatest_elem];
  b[row_with_greatest_elem] = temp;


}

/**
 * \brief forward_elim
 * \decription transforms the matrix A into the _row echelon form_ wherein the
 * (i,j)th term is reduced to 0,  and all entries in a column below the leading
 * coefficient are zeros.
 * \params n, the size of the matrix
 *         A, the LHS of the matrix, filled with random variables
 *         b, the RHS of the matrix, here as a vector of answers
 * \return void function
 */
void forward_elim(int n, double* A, double* b ){
  //want to solve for the x[n-1] variable and the A[n-1] row

  //using the ith row to eliminate the ith item from the jth row
  for( int i = 0; i < n; i++){
    //pivot(n, A, b, i, thread_count);
    pivot(n, A, b, i);
#   pragma omp parallel for 
    for (int j = i+1; j < n; j++){
      //      if(omp_get_thread_num() == 0){
      //        printf("number of threads = %d\n", omp_get_num_threads());
      //            }
      double scaling_factor = (A[j*n+i] / A[i*n+i]);
      for(int k = i; k < n; k++){
        //everything to the left of the ith item is eliminated already
        A[j*n+k] -= A[i*n+k] * scaling_factor;
      }
      b[j] -= b[i] * scaling_factor;
    }
    //gauss
    //use row[i] to eliminate variable[1][2] in all rows > i
  }
}

/**
 * \brief back substituion
 * \decription we can assume at this point that the matrix we are dealing with
 * has been swapped/pivoted correctly and is now in the Upper Triangular form,
 * meaning that the matrix would have the variables in this form:
 *                    [ x  x  x ]
 *                    [ 0  x  x ]
 *                    [ 0  0  x ]
 *  once it is that form, we can solve the x_3 value  and then work backward to
 *  solve for all the unknowns
 * \params n, the size of the matrix
 *         A, the LHS of the matrix, filled with random variables
 *         b, the RHS of the matrix, here as a vector of answers
 *         x, a vector of unknowns
 * \return void function
 */
void back_substitution(int n, double* A, double* b, double* x){

  //backward loop
  for (int i = n-1; i >= 0; i--){
    //solve from bottom right corner of matrix
    double new_b = b[i];
    for (int j = n-1; j >= i; j--){
      // j = i === diagonal, which will have nonzero values
      // shifting the values of A that are not on the diagonal over to b to
      // solve for x
      if (i == j){
        x[i] = new_b / A[i*n+j];
      }
      new_b -= A[i*n+j] * x[j];
    }
  }
}

/**
 * \brief main driver
 * \description initializes variables, matrix, vectors. feeds to forward_elim
 * and then back_substituion 
 * \params command line arguments
 * \return 0
 */
int main(int argc, char ** argv){

  int n;
  if (argc != 2){
    printf("./gepp <size of n>\n ");
    exit(1);
  }


  n = atoi(argv[1]);


  double start, end;

  start = omp_get_wtime();

  srand48(time(0));

  int num_procs = omp_get_num_procs();
  printf("num procs = %d\n", num_procs);

  double *A;
  double *A_orig;
  double *b;
  double* b_orig;
  double *x;

  A = (double*)malloc(n*n*sizeof(double)); // matrix of doubles of size n
  A_orig = (double*) malloc(n*n*sizeof(double));
  b =  (double*) malloc(n*sizeof(double));
  b_orig =  (double*) malloc(n*sizeof(double)); //for Lnorm and printing
  x = (double*)malloc(n* sizeof(double));// unknown variables

  if (n <11){
    // read the augmented matrix for A | b from stdio, one row per line
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n+1; j++){

        double number;
        scanf("%lf", &number);
        if (j < n){
          A[i*n+j] = number;
          A_orig[i*n+j] = number;
        }
        b[i] = number;
        b_orig[i] = number;
      }
    }// A array values
  }else{
    /**fill A | b with random numbers using drand48 random num generator
     * scale result to the range of [-1.0e6, 1.0e6) if on the main diagonal
     * scale result to [-1.03e, 1.03e) if not on the diagonal.
     */
    for(int i = 0; i < n; i ++){
      for (int j = 0; j < n; j++){
        double number = drand48(); // gets a number [0,1)
        number -= 0.5; // gets a numnber [-.5, .5)
        number = number * 2000000; // gets a number [-1.0e6, 1.0e6)

        A[i*n+j] = number;
        A_orig[i*n+j] = number;
      }
    } // filling A array with random values

    for(int i = 0; i < n; i++){
      double number = drand48();
      number -= 0.5;
      number = number * 2000000;
      b[i] = number;
      b_orig[i] = number;
    } // filling b array with random values
  }

  forward_elim(n, A, b);
  back_substitution(n, A, b, x);


  double norm = L2norm(n, A_orig, b_orig, x);

  end = omp_get_wtime();

  printf("elapsed time = %f seconds\n",  end - start);

  if (n < 11){
    for (int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        printf("%0.10e ", A_orig[i*n+j]);
      }
      printf("%0.10e\n", b_orig[i]);
    }

    for(int i = 0; i < n; i++){
      printf("%0.10e ", x[i]);
    }
    printf("\n");
  }

  printf("L2-norm = %0.10e\n", norm);

  free(A);
  free(A_orig);
  free(b);
  free(b_orig);
  free(x);

  return 0;
}

