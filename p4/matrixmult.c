// p4
// Matrix Multiplication with IJK forms and MPI
// Jeff Anderson
//
// with resources:
// https://stackoverflow.com/questions/27575912/mpi-debugging-segmentation-fault
// https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node70.html
// a good example of how to use MPI_Scatterv:
//          https://gist.github.com/ehamberg/1263868/cae1d85dee821d45fb0cd58747aaf33370f3f1ed
// and code provided by Dr. Challinger and the textbook
// An Introduction to Parallel Programming by Peter Pacheco
//


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h> // for string comparison
#include <mpi.h>  // for MPI functionality

/**
 * @func: ijk
 * @desc: performs matrix multiplication in typical ijk form
 * @return: void func
 */
void ijk(int n, int* A, int* B, int* C, int sendcnts)
{
  int i, j, k;
  int size = sendcnts / n;

  for (i = 0; i < size; i++)
  {
    for (j = 0; j < n; j++)
    {
      for (k = 0; k < n; k++)
      {
        //C[i][j] = C[i][j] + A[i][k] * B[k][j];
        C[i*n + j] = C[i*n + j] + A[i*n + k] * B[k*n + j];
      }
    }
  }
}

/**
 * @func: ikj
 * @desc: performs matrix multiplication in ikj form
 * @return: void func
 */
void ikj(int n, int* A, int* B, int* C, int sendcnts)
{

  int i, j, k;
  /* while typical matrix multiplication is in [i , j , k] form, that is 
   * not always the best way to vectorize efficiently.
   * Interchanging the loops make all the references STRIDE-1
   * as referenced in the Intel compiler docs to make the program more effcient
   */
  int size = sendcnts / n;

  for (i = 0; i < size; i++)
  {
    for (k = 0; k < n; k++)
    {
      for (j = 0; j < n; j++)
      {
        C[i * n + j] = C[i * n + j] + A[i * n + k] * B[k * n + j];
      }
    }
  }
}

/**
 * @func: kij
 * @desc: multiplies two matrices in kij form
 * @return: void func
 */
void kij(int n, int* A, int* B, int* C, int sendcnts)
{
  int i, j, k;
  int size = sendcnts / n;
  for(int x = 0; x < sendcnts; x++){
    C[x] = 0;
  }

  for (k = 0; k < n; k++)
  {
    for (i = 0; i < size; i++)
    {
      for (j = 0; j < n; j++)
      {
        C[i*n + j] = C[i*n + j] + A[i*n + k] * B[k*n + j];

      }
    }
  }
}

/*
 * @func my_rand
 *        a random number generator between 100, 0
 * @return returns a random number between 0-99
 */
int my_rand()
{
  // return a random int between 100 and 0
  return rand() % 101;
}

/**
 * @func main
 * @desc: main driver of the program
 *    takes in 2 user-inputted or randomized matrices and produces the dot
 *    product using MPI function calls and parallelization.
 * @return: nothing
 */
int main()
{

  int comm_sz; // number of processes
  int my_rank; // variable to tell which process is which

  char form[4]; // 4 chars to leave space for null terminator
  char flag;
  int n;
  double start_time, finish_time; // these have to be a double per the man pages

  int *A; // set up A, B, C
  int *B;
  int *C;

  MPI_Init(NULL, NULL); // initialize MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); // grab number of processes/cores
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // find out how many processes are being used

  // set up in rank 0 for all the other cores to us
  if (my_rank == 0){
    printf("enter form, then flag, then n \n");
    scanf("%3s\n", form); // read in only 3 characters - scanf will add a null terminator
    scanf("%c", &flag); // read in [r]andom or [i]nput
    scanf("%d", &n); // read in size of matrix

    printf("running on %d processor(s)\n", comm_sz);

    //initialze A, B, C and set to 0
    A = (int*)calloc(n*n, sizeof(int));
    B = (int*)calloc(n*n, sizeof(int));
    C = (int*)calloc(n*n, sizeof(int));


    // for determining user choice
    if (flag == 'R')
    {

      // random generation of matrices
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          A[i*n + j] = my_rand();
        }
      } // A for loop

      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          B[i * n + j] = my_rand();
        }
      } // B for loop
    }
    else if (flag == 'I')
    {
      // input matrices
      // set input values for A array
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          int number;
          scanf("%d", &number);
          A[i * n + j] = number;
        }
      } // A array values

      // set input values for B array
      for (int i = 0; i < n; i++)
      {
        for (int j = 0; j < n; j++)
        {
          int number;
          scanf("%d", &number);
          B[i * n + j] = number;
        }
      } // B array values
    }
    else
    {
      //error check
      fprintf(stderr, "please choose either [I]nput matrices or choose [R]andom matrices\n");
      exit(EXIT_FAILURE); // stdlib defines EXIT_FAILURE as 1, but is clearer in code
    }
  }


  MPI_Barrier(MPI_COMM_WORLD); 
  // after all the initial set up, this is the first MPI call.
  // per the man page:
  // Blocks the caller until all processes in the communicator have called it;
  // that is, the call returns at any process only after all members of the
  // communicator have entered the call.
  if (my_rank ==  0){
    // set up stopwatch
    start_time = MPI_Wtime();
  }

  // bcast sizeof matrices to all processes
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int *sendcnts = (int*)calloc(comm_sz, sizeof(int)); // how many elements to send to each process
  // using calloc instead of malloc because calloc will allocate the size we
  // need and also initialize all the values to 0 at the same time. 
  int *displacements = (int*)calloc(comm_sz, sizeof(int));
  int remainder = n % comm_sz; // elements remaining after division among processes
  int sum = 0; // used for calculating displacements -- sum of counts


  // calculate send counts and displacements
  for (int i = 0; i < comm_sz; i++){
    sendcnts[i] = (n / comm_sz) * n;
    if (remainder > 0){
      remainder--;
      sendcnts[i]++;
    }
    displacements[i] = sum;
    sum += sendcnts[i];
  }

  //initialize a local a & c for each individual process to use
  int *local_a = (int *) calloc(sendcnts[my_rank], sizeof(int));
  int *local_c = (int *) calloc(sendcnts[my_rank], sizeof(int));

  if(my_rank != 0){
    A = (int*)calloc(n*n, sizeof(int));
    B = (int*)calloc(n*n, sizeof(int));
    C = (int*)calloc(n*n, sizeof(int));
  }
  // bcast which form to all processes
  MPI_Bcast(&form, 4, MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast(B, n * n, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatterv(A, sendcnts, displacements, MPI_INT, local_a, sendcnts[my_rank],
      MPI_INT, 0, MPI_COMM_WORLD);

  // for determining form
  if (strcmp(form, "ijk") == 0)
  {
    ijk(n, local_a, B, local_c, sendcnts[my_rank]);
  }
  else if (strcmp(form, "ikj") == 0)
  {
    ikj(n, local_a , B, local_c, sendcnts[my_rank]);
  }
  else if (strcmp(form, "kij") == 0)
  {
    kij(n, local_a, B, local_c, sendcnts[my_rank]);
  }
  else
  {
    //error check
    fprintf(stderr, "please choose ijk, ikj, or kij");
    exit(EXIT_FAILURE); // stdlib defines EXIT_FAILURE as 1, but is clearer in code
  }

  // Gatherv works in conjucntion with Scatterv 
  MPI_Gatherv(local_c, sendcnts[my_rank], MPI_INT, C, sendcnts, displacements, MPI_INT, 0, MPI_COMM_WORLD);

  if (my_rank == 0){
    printf("form = %s\n", form); // print the form
    finish_time = MPI_Wtime(); // stop stopwatch
    printf("elapsed time = %f seconds\n", finish_time - start_time);
  }

  // print calculated dot product if inputted
  if (flag == 'I')
  {

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        //printf("%d ", C[i][j]);
        printf("%d ", C[i*n+j]);
      }
      printf("\n");
    }
  }

  //free all allocated memory
  free(A);
  free(B);
  free(C);
  free(local_a);
  free(local_c);
  free(sendcnts);
  free(displacements);

  //end MPI
  MPI_Finalize();

  return 0;
}
