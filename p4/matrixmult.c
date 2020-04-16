// p4

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h> // for string comparison
//#include <mpi.h>  // for MPI functionality

void ikj(int n, int A[n][n], int B[n][n], int C[n][n])
{

    int i, j, k;
    /* while typical matrix multiplication is in [i , j , k] form, that is 
   * not always the best way to vectorize efficiently.
   * Interchanging the loops make all the references STRIDE-1
   * as referenced in the Intel compiler docs to make the program more effcient
   */
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < n; k++)
        {
            for (j = 0; j < n; j++)
            {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

void ijk(int n, int A[n][n], int B[n][n], int C[n][n])
{
    int i, j, k;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < n; k++)
            {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

void kij(int n, int A[n][n], int B[n][n], int C[n][n])
{
    int i, j, k;

    for (k = 0; k < n; k++)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

/*
 * @func my_rand
 *        a random number generator between 100, 0
 * @return returns a random number 
 */
int my_rand()
{
    // return a random int between 100 and 0
    return rand() % 101;
}

int main()
{

    int comm_sz; // number of processes
    int my_rank; // variable to tell which process is which

    char form[4]; // 4 chars to leave space for null terminator
    char flag;
    int n;

    printf("enter form, then flag, then n \n");
    scanf("%3s\n", form); // read in only 3 characters - scanf will add a null terminator
    scanf("%c", &flag);
    scanf("%d", &n);

    int A[n][n];
    int B[n][n];
    int C[n][n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[i][j] = 0;
        }
    }

    srand(time(0));

    //MPI_Init
    // for determining user choice
    if (flag == 'R')
    {

        // random generation of matrices
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = my_rand();
            }
        } // A for loop

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                B[i][j] = my_rand();
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
                A[i][j] = number;
            }
        } // A array values

        // set input values for B array

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int number;
                scanf("%d", &number);
                B[i][j] = number;
            }
        } // B array values
    }
    else
    {
        //error check
        fprintf(stderr, "please choose either [I]nput matrices or choose [R]andom matrices");
        exit(EXIT_FAILURE); // stdlib defines EXIT_FAILURE as 1, but is clearer in code
    }
    printf("form = %s\n", form);
    // for determining form
    if (strcmp(form, "ijk") == 0)
    {
        ijk(n, A, B, C);
    }
    else if (strcmp(form, "ikj") == 0)
    {
        ikj(n, A, B, C);
    }
    else if (strcmp(form, "kij") == 0)
    {
        kij(n, A, B, C);
    }
    else
    {
        //error check
        fprintf(stderr, "please choose ijk, ikj, or kij");
        exit(EXIT_FAILURE); // stdlib defines EXIT_FAILURE as 1, but is clearer in code
    }

    printf("running on %d processor\n", 1);
    printf("elapsed time = %d seconds\n", 0);

    if (flag == 'I')
    {

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                printf("%d ", C[i][j]);
            }
            printf("\n");
        }
    }
}
