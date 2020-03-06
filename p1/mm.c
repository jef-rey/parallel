/*
 * 551 - program 1 - AUTOMATIC VECTORIZATION
 * Jeff Anderson 008326891
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h> // for getruage and gettimeofday
#include <sys/time.h> // for getruage and gettimeofday
#include <time.h> // for srand()


/*
 * @fucn add
 * adds two matrices of nxn size 
 */
void add(int n, float A[n][n], float B[n][n], float C[n][n]){
  // avoiding int **A as per intel c++ compiler guidelines

  int i, j;
  for (i = 0; i < n; i++){
#pragma vector aligned
    for (j = 0; j < n; j++){
      C[i][j] =  A[i][j] + B[i][j];
    }
  }
}

/*
 * @func multiply
 * multiplies two matrices of nx size using ikj (aka stride-1)
 * for faster optimization
 */
void multiply(int n, float A[n][n], float B[n][n], float C[n][n]){

  int i, j, k;
  /* while typical matrix multiplication is in [i , j , k] form, that is 
   * not always the best way to vectorize efficiently.
   * Interchanging the loops make all the references STRIDE-1
   * as referenced in the Intel compiler docs to make the program more effcient
   */
  for (i = 0; i < n; i++){
#pragma vector aligned 
    for (k = 0; k < n; k++){
      for (j = 0; j < n; j++){
        C[i][j] = C[i][j] + A[i][k] * B[k][j];
      }
    }
  }
}


/*
 * @func my_rand
 *        a random number generator between -50, 50
 * @return returns a random number 
 */
float my_rand(){

  float random = (float)rand() / RAND_MAX;
  // this generates a random number between 0 -1
  random = random * 100 - 50;
  // then multiply by 100 and subtract 50 to be within bounds [-50,50]


  return random;
}


/*
 * @func main
 *        driver for the program
 */
int main(int argc, char *argv[]){

  struct rusage usage; // for getrusage()
  struct timeval start_add; // for gettimeofday()
  struct timeval start_mult;
  struct timeval end_add; // clearer to me in code to have a start/end for add/mult
  struct timeval end_mult;

  int n;
  char choice;
  //char* choice[2]; // size 2- one for char, one for space. strings in C
  // are null-terminated character arrays.
  //scanf("%c", choice);
  scanf("%c", &choice);
  scanf("%d", &n);

  srand(time(0));

  float random =  (rand() % (50 - (-50) + 1)) - 50;


  /* 
   * telling the compiler to align A, B, & C arrays on 16-bye boundaries
   * this allows us to use aligned instructions and produce faster code.
   */
#ifdef ALIGNED
  float A[n][n]   __attribute__((aligned(16))); // declaring A
  float B[n][n]   __attribute__((aligned(16)));
  float C_mult[n][n]  __attribute__((aligned(16)));
  // start C off initialized to 0. This array 
  // will be used for multiplication
  float C_add[n][n]   __attribute__((aligned(16))); // this C will be used for addition.
#endif

#pragma novector
  for(int q = 0; q < n; q++){
    for(int z = 0; z < n; z++){
      C_mult[q][z] = 0;
      C_add[q][z] = 0;
    }
  }


  if (choice == 'R'){
    //randomize matrices

#pragma novector
    for (int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        A[i][j] = my_rand();
      }
    }// A for loop


#pragma novector
    for (int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        B[i][j] = my_rand();
      }
    }// B for loop

  }else if (choice == 'I'){
    // input matrices
    // set input values for A array
#pragma novector
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        float number;
        scanf("%f", &number);
        A[i][j] = number;
      }
    }// A array values

    // set input values for B array
#pragma novector
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        float number;
        scanf("%f", &number);
        B[i][j] = number;
      }
    }// B array values


  }else{
    fprintf(stderr, "please choose either [I]nput matrices or choose [R]andom matrices");
    exit(EXIT_FAILURE); // stdlib defines EXIT_FAILURE as 1, but is clearer in code
  }

  gettimeofday(&start_add, NULL);
  add(n, A, B, C_add);
  gettimeofday(&end_add, NULL);


  gettimeofday(&start_mult, NULL);
  multiply(n, A, B, C_mult);
  gettimeofday(&end_mult, NULL);

  printf("\n");


  if (choice == 'I'){
#pragma novector
    for(int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        printf("%.0f ", C_add[i][j]);
      }
      printf("\n");
    }

    printf("\n");
#pragma novector
    for(int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        printf("%.0f ", C_mult[i][j]);
      }
      printf("\n");
    }
  }



  getrusage(RUSAGE_SELF, &usage);
  unsigned long user = usage.ru_utime.tv_sec * 1000000 + usage.ru_utime.tv_usec; 

  float add_time = (end_add.tv_sec + end_add.tv_usec /1000000.0)-( start_add.tv_sec + start_add.tv_usec /1000000.0) ;
  float mult_time = (end_mult.tv_sec + end_mult.tv_usec /1000000.0) -(start_mult.tv_sec + start_mult.tv_usec / 1000000.0);

  printf("time for adding matrices = %f\n", add_time);
  printf("time for multiplying matrices = %f\n", mult_time);

  printf("usage time = %ld", user);
  printf("\n");
  printf("Max RSS = %ld", usage.ru_maxrss);
  printf("\n");

  exit(EXIT_SUCCESS); //exit_success === exit(0) but clearer to read in code
}
