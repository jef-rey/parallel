// p5
// Gaussian Elimination with Partial Pivoting
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

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

void pivot(int n, double* A, double* b, int diag ){
  int i;
  
  double greatest_val = 0;
  int row_with_greatest_elem = diag;
  
/*#pragma omp parallel for num_threads(thread_count) \
  private(i) shared(n, A, b, diag) */
  for (/*int*/ i = diag; i < n; i++){
    double abs_elem = fabs(A[i*n+diag]);
    if (abs_elem > greatest_val){
      greatest_val = abs_elem;
      row_with_greatest_elem = i;
    }
  }

  //swap rows in A column
/*#pragma omp parallel for num_threads(thread_count) \
  private(i) shared(n, A, b, diag)*/
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

void back_substitution(int n, double* A, double* b, double* x){

  //backward loop
  for (int i = n-1; i >= 0; i--){
    //solve from bottom right corner of matrix
    double new_b = b[i];
    for (int j = n-1; j >= i; j--){
      // j = i === diagonal, which will have nonzero values
      // shifting the values of A that are not on the diagonal over to B to
      // solve for x
      if (i == j){
        x[i] = new_b / A[i*n+j];
      }
      new_b -= A[i*n+j] * x[j];
    }
  }
}

int main(int argc, char ** argv){

  int n;
  if (argc != 2){
    printf("./gepp <size of n>\n ");
    // this is more a relic of testing, since we ran this using the sbatch
    // system
    exit(1);
  }


  n = atoi(argv[1]);
  

  double start, end;

  start = omp_get_wtime();

  srand48(time(0));

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

