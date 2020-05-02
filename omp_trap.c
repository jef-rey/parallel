/**
 * program 3 -- REVISED -- FROM MPI TO OMP
 * parallel numical integration with open mpi
 * Jeff Anderson
 * April 2020
 *
 * to compile:
 *  gcc -o omp_trap -fopenmp omp_trap.c -lm
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h> // for cosl, sinl, fabsl, and powl
#include <omp.h>  // for MPI functionality  

/**
 * @func: f
 * @desc: the integral provided for the basis of this project
 * @return: returns f(x) for a provided x
 */
long double f(long double x){
  return (-6*cosl(x/9) +sqrtl(powl((5*sinl(x/3) + 8*sinl(x/2)), 4)) + 10);
}

/**
 * @func: approx_area
 * @desc: takes user input for lower bound, upper bound, and number of
 * trapezoids 
 * @return: returns the sum of the areas of the trapezoids
 *
 * adapted from "Introduction of Parallel Programming" by Peter Pacheco
 */
long double approx_area(long double a, long double b, int n){

  long double h = (b - a) / n; // determine trap width
  // n is the number of trapezoids we're using at the particular iteration

  long double local_a; // going to use local_a as the next point to iterate on.

  long double sum = (f(a) + f(b)) / 2;
  // adapted from pacheco book
#  pragma omp parallel for num_threads(2) \
  reduction(+: sum)
  for (int i = 1; i <= n-1; i++){
    local_a = a + i * h;
    sum += f(local_a);
  }

  sum *= h;

  return sum;
}

/**
 * @func: print_report
 * @desc: prints the output as outlined in the documentaion for the assignment
 * @return: void
 */
void print_report(double start_time, double finish_time, int n, long double a, long double b, long double sum ){

  long double stop_criteria = 0.5E-14;
  long double true_value = 6745.9948824470072783;	//	true value of integral (using 20 digits)
  long double arte = fabsl((true_value - sum) / true_value); // absolute relative true error
  /*
   * ARTE = abs( relative true error)
   * RTE = (true error / true value)
   */

  printf("Elapsed time = %.6e seconds\n", finish_time - start_time);
  printf("With n = %.0d trapezoids, our estimate\n", n);
  printf("of the integral from %Lf to %Lf = %.13Le\n", a, b, sum);
  printf("true value = %.19Le\n", true_value);
  printf("absolute relative true error = %.19Le\n", arte);
  if (arte < stop_criteria)
    printf(" is less than criteria = %.19Le\n", stop_criteria);
  else if (arte == stop_criteria)
    printf(" is equal to the criteria = %.19Le\n", stop_criteria);
  else
    printf(" is NOT less than criteria = %.19Le\n", stop_criteria);
}


/**
 * @func: main
 * @desc: driver for program. 
 * Initializes the Message Paasing Interface (MPI) to evenly distribute the work
 * amongst multiple cores on multiple lab machines.
 *
 * files mpi_trap2.c , mpi_trap3.c , and mpi_trap4.c were heavily source for
 * figuring out how to implement MPI processess correctly. 
 * as was :
 * https://computing.llnl.gov/tutorials/mpi/
 * 
 */
int main(void){

  int comm_sz; // number of processes -- same as paecheo used. 
  int my_rank;  // my rank

  long double a, b;	// a = lower bound, b = upper bound
  int n; // number of trapezoids

  long double h;	//	width of each trapezoid
  long double sum;	//	Global sum of all trapezoid areas
  double start_time, finish_time;

  //MPI_Init(NULL, NULL); // initialize MPI

  //MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); // grab number of processes/cores

  /* Find out how many processes are being used */
  //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // grab user input
  if (my_rank == 0)
  {
    printf("Enter a, b and n\n");
    scanf("%Lf %Lf %d", &a, &b, &n);
    printf("\nRunning on %d cores.\n", comm_sz);
  }

  //MPI_Barrier(MPI_COMM_WORLD);// MPI_BARRIER creates a barrier synchronization in a group. 
  // each task blocks until all the tasks in the group reaches the same call,
  // then all tasks proceed.
  if (my_rank == 0)
  {
    //start_time = MPI_Wtime(); // start timer
    start_time = omp_get_wtime(); // start timer
  }


  h = (b - a) /  n; // get width of trapezoids

  long double local_a, local_b, local_sum;
  int local_n;

  // length of each processes interval of integration
  local_n = n / comm_sz;				// number of trapezoids for the  process
  local_a = a + my_rank*local_n*h; 	// lower bound for the processes
  local_b = local_a + local_n*h;		// upper bound for the processes 


  local_sum = approx_area(local_a, local_b, local_n);

  // add all the sums and send to first process -- like that tree graph
  //MPI_Reduce(&local_sum, &sum, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (my_rank == 0)
  {
    //finish_time = MPI_Wtime(); // stop timing
    finish_time = omp_get_wtime(); // stop timing

    print_report(start_time, finish_time, n, a, b, sum); // calls print_report to print report
  }

  //MPI_Finalize(); // stop MPI communication
  return 0;
}
