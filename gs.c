#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* By Jason Yao */
typedef enum { false, true } bool; // Enables boolean types

/* Configuration settings */
bool IS_DEBUG_MODE = false;         /* Change to true to see intermediate values */
bool IS_SEQUENTIAL_MODE = false;     /* Change to true to switch to sequential version */

/***** Globals ******/
float **a;                          /* The coefficients */
float *x;                           /* The unknowns */
float *b;                           /* The constants */
float *curr;                        /* current values of unknown */
float err;                          /* The absolute relative error */
int num = 0;                        /* number of unknowns */
bool SOLUTION_IS_SOLVED = false;    /* Whether the solution is already solved */

/* Original function declarations */
void check_matrix();    /* Check whether the matrix will converge */
void get_input();       /* Read input from file */

/* Solution suite function declarations */
void generateNextRound();
int sequential();
void cleanup();

/* Testing suite function declarations */
void testCurrent();
void testAll();
void testGeneratedXAndErr(float *newX_i, float *newErr);

/**********************
 * Test Suite Functions
 **********************/
void testCurrent()
{
    printf("The curr[i] matrix:\n");
    for (int i = 0; i < num; ++i)
    {
        printf("%f ", curr[i]);
    }
    printf("\n");
    printf("---------------------------------------------\n");
} // End of test current function

void testAll()
{
    printf("The a[i][i] matrix:\n");
    for (int i = 0; i < num; ++i)
    {
        for (int j = 0; j < num; ++j)
            printf("%f ", a[i][j]);
        printf("\n");
    }
    printf("---------------------------------------------\n");

    printf("The x[i] matrix:\n");
    for (int i = 0; i < num; ++i)
    {
        printf("%f ", x[i]);
    }
    printf("\n");
    printf("---------------------------------------------\n");

    testCurrent();

    printf("The b[i] matrix\n");
    for (int i = 0; i < 3; ++i)
        printf("%f ", b[i]);
    printf("\n");
    printf("---------------------------------------------\n");

    printf("Error: %f\n\n", err);
} // End of the test all function

void testGeneratedXAndErr(float *newX_i, float *newErr)
{
    printf("-------------------------\n");
    for (int i = 0; i < num; ++i)
        printf("Newly generated x_i: %f\n", newX_i[i]);
    printf("-------------------------\n");

    printf("-------------------------\n");
    for (int i = 0; i < num; ++i)
        printf("Newly generated Errors: %f\n", newErr[i]);
    printf("-------------------------\n");
} // End of the testGeneratedXAndErr function

/**************************
 * Solution Suite Functions
 **************************/
void generateNextRound()
{
    float* newX_i = (float *) malloc(num * sizeof(float));
    float* newErr = (float *) malloc(num * sizeof(float));

    // Generates all new x_i
    for (int i = 0; i < num; ++i)
    {
        // Finds x_i
        float numerator = b[i];
        for (int j = 0; j < num; ++j)
            if (j != i)
                numerator -= (a[i][j] * curr[j]);

        float denominator = a[i][i];
        float x_i = numerator / denominator;

        newX_i[i] = x_i;
    }

    // Checks each x_i for closeness validity
    bool isDone = true;
    for (int i = 0; i < num; ++i)
    {
        newErr[i] = fabs((newX_i[i] - curr[i])/(newX_i[i]));
        if (newErr[i] > err)
            isDone = false;

        curr[i] = newX_i[i];
        x[i] = curr[i];
    }

    if (IS_DEBUG_MODE)
        testGeneratedXAndErr(newX_i, newErr);

    if (isDone)
        SOLUTION_IS_SOLVED = true;

    free(newX_i);
    free(newErr);
} // End of the generate next round function

/**
 * Parallel version of the code to solve n equations to within a bounds
 */
int parallel()
{
    // Sets the current x values to the initial x values
    for (int i = 0; i < num; ++i)
        curr[i] = x[i];

    /* MPI stuff */
    int comm_size;      // Number of processes
    int my_rank;        // My process rank
    int count = 0;      // The current iteration cycle

    /**
     * Solution:
     * 1.) send to each process n/p_ceiling old x_i values in a subarray,
     * where n == number of elements of x_i, and p == number of processes.
     * 2.) Have each process calculate the new x_i subarray
     *
     * 3.) Use MPI_Allgather to gather all new x_i subarrays. For each new x_i, calculate new_error.
     * For each element in the full error subarray, check to see if any are above the absolute error stated.
     * If no, sets x to the new x_i values, and returns the count. If yes, starts from 1 with the MPI_Allgather.
     */

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int local_n = (num + comm_size - 1)/(comm_size); // sub_array_size_ceiling, number of elements in the subarray
    float* local_curr; // recieving buffer for elements being scattered in blocks
    float* new_x = (float *) malloc(num * sizeof(float));       // Contains all new values
    float* new_err = (float *) malloc(num * sizeof(float));     // Contains all new errors

    for (; !SOLUTION_IS_SOLVED; ++count)
    {
        local_curr = (float *) malloc(local_n * sizeof(float)); // Buffer for calculated x_i subarray

        if (my_rank == 0)
        {
            // Process is root
            // Scatters the current x_i values to each process
            MPI_Scatter(curr, local_n, MPI_FLOAT, local_curr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        } // End of dealing with root rank
        else
        {
            // Process is not root
            MPI_Scatter(a, local_n, MPI_FLOAT, local_curr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
        } // End of dealing with non-root rank

        // Generates all new x_i
        for (int i = 0; i < local_n; ++i)
        {
            // Finds x_i
            local_curr[i] = b[i];

            for (int j = 0; j < num; ++j)
                if (j != i)
                    local_curr[i] -= (a[i * num][j] * curr[j]);
            local_curr[i] /= a[i][i];
        } // At the end of this, local_x[] holds the newly calculated x_i subarray

        // Gathers each processes's local_x[], concatenates, and stores in each process's curr[]
        MPI_Allgather(local_curr, local_n, MPI_FLOAT, new_x, local_n, MPI_FLOAT, MPI_COMM_WORLD);

        // Checks each x_i for closeness validity
        bool isDone = true;
        for (int i = 0; i < num; ++i)
        {
            new_err[i] = fabs((new_x[i] - curr[i])/(new_x[i]));
            if (new_err[i] > err)
                isDone = false;

            curr[i] = new_x[i];
            x[i] = curr[i];
        }

        if (IS_DEBUG_MODE)
            testGeneratedXAndErr(new_x, new_err);

        if (isDone)
            SOLUTION_IS_SOLVED = true;
    } // End of parallel

    // Frees everything
    free(local_curr);
    free(new_x);
    free(new_err);




    // Gathers all of the newly calculated x_i values
//    MPI_Allgather(local_x, sub_array_size_ceiling, MPI_FLOAT, curr, sub_array_size_ceiling, MPI_FLOAT, MPI_COMM_WORLD);

//    // Version based off of allgather
//    local_x = (float *) malloc(num * sizeof(float));    // Out
//    new_err = (float *) malloc(num * sizeof(float));    // Out
//    int local_m, local_n; // m == upper range for i, n == upper range for j
//    int local_i, local_j;
//
//    for (; !SOLUTION_IS_SOLVED; ++count)
//    {
//        // Generates in a parallel fashion
//        MPI_Allgather(local_x, local_n, MPI_FLOAT, curr, local_n, MPI_FLOAT, MPI_COMM_WORLD);
//
//        // Generates all new x_i
//        for (local_i = 0; local_i < local_m; ++local_i)
//        {
//            // Finds x_i
//            local_x[local_i] = b[local_i];
//
//            for (local_j = 0; local_j < num; ++local_j)
//                if (local_j != local_i)
//                    local_x[local_i] -= (a[local_i*num][local_j] * curr[local_j]);
//
//            local_x[local_i] /= a[local_i][local_i];
//        }
//
//        // Checks each x_i for closeness validity
//        bool isDone = true;
//        for (int i = 0; i < num; ++i)
//        {
//            new_err[i] = fabs((local_x[i] - curr[i])/(local_x[i]));
//            if (new_err[i] > err)
//                isDone = false;
//
//            curr[i] = local_x[i];
//            x[i] = curr[i];
//        }
//
//        if (IS_DEBUG_MODE)
//            testGeneratedXAndErr(local_x, new_err);
//
//        if (isDone)
//            SOLUTION_IS_SOLVED = true;
//        /* End of parallel code */
//
//        if (SOLUTION_IS_SOLVED)
//            break;
//    } // End of iterating through all test cycles





    // Version based off of rank differential
//    if (my_rank != 0)
//    {
//        new_x_i = (float *) malloc(num * sizeof(float));
//        new_err = (float *) malloc(num * sizeof(float));
//
//        // Process is not root
//        for (; !SOLUTION_IS_SOLVED; ++count)
//        {
//
//            generateNextRound();
//
//            if (SOLUTION_IS_SOLVED)
//                break;
//        }
//
//    } // End of dealing with root rank
//    else
//    {
//        // Process is root
//    } // End of dealing with non-root rank

//    free(local_x);
//    free(new_err);
    MPI_Finalize();
    return count;
} // End of the parallel function

/**
 * Returns the number of iterations to solve these equations
 */
int sequential()
{
    // Sets the current x values to the initial x values
    for (int i = 0; i < num; ++i)
        curr[i] = x[i];

    int count = 0;
    for (; !SOLUTION_IS_SOLVED; ++count)
    {
        generateNextRound();

        if (SOLUTION_IS_SOLVED)
            break;
    }
    return count;
} // End of the sequential function

/* 
   Conditions for convergence (diagonal dominance):
   1. diagonal element >= sum of all other elements of the row
   2. At least one diagonal element > sum of all other elements of the row
 */
void check_matrix()
{
    int bigger = 0; /* Set to 1 if at least one diag element > sum  */
    int i, j;
    float sum = 0;
    float aii = 0;

    for(i = 0; i < num; i++)
    {
        sum = 0;
        aii = fabs(a[i][i]);

        for(j = 0; j < num; j++)
            if( j != i)
                sum += fabs(a[i][j]);

        if( aii < sum)
        {
            printf("The matrix will not converge\n");
            exit(1);
        }

        if(aii > sum)
            bigger++;

    }

    if( !bigger )
    {
        printf("The matrix will not converge\n");
        exit(1);
    }
} // End of the check matrix function

/******************************************************/
/* Read input from file */
void get_input(char filename[])
{
    FILE * fp;
    int i,j;

    fp = fopen(filename, "r");
    if(!fp)
    {
        printf("Cannot open file %s\n", filename);
        exit(1);
    }

    fscanf(fp,"%d ",&num);
    fscanf(fp,"%f ",&err);

    /* Now, time to allocate the matrices and vectors */
    a = (float**)malloc(num * sizeof(float*));
    if( !a)
    {
        printf("Cannot allocate a!\n");
        exit(1);
    }

    for(i = 0; i < num; i++)
    {
        a[i] = (float *)malloc(num * sizeof(float));
        if( !a[i])
        {
            printf("Cannot allocate a[%d]!\n",i);
            exit(1);
        }
    }

    x = (float *) malloc(num * sizeof(float));
    if( !x)
    {
        printf("Cannot allocate x!\n");
        exit(1);
    }

    curr = (float *) malloc(num * sizeof(float));
    if( !curr)
    {
        printf("Cannot allocate curr!\n");
        exit(1);
    }

    b = (float *) malloc(num * sizeof(float));
    if( !b)
    {
        printf("Cannot allocate b!\n");
        exit(1);
    }

    /* Reading in from file */

    /* The initial values of Xs */
    for(i = 0; i < num; i++)
        fscanf(fp,"%f ", &x[i]);

    for(i = 0; i < num; i++)
    {
        for(j = 0; j < num; j++)
            fscanf(fp,"%f ",&a[i][j]);

        /* reading the b element */
        fscanf(fp,"%f ",&b[i]);
    }
    fclose(fp);
} // End of the get input function


/************************************************************/

void cleanup()
{
    for (int i = 0; i < num; ++i)
        free(a[i]);
    free(a);
    free(x);
    free(curr);
    free(b);
} // End of the cleanup function

int main(int argc, char *argv[])
{
    /* Timing stuff */
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    /* Number of iterations */
    int nit = 0;

    if( argc != 2)
    {
        printf("Usage: ./gsref filename\n");
        exit(1);
    }

    /* Read the input file and fill the global data structure above */
    get_input(argv[1]);

    /* Check for convergence condition */
    check_matrix();

    /* Solves the matrix */
    if (IS_SEQUENTIAL_MODE)
        nit = sequential();
    else
        nit = parallel();

    /* Writing to the stdout */
    for(int i = 0; i < num; i++)
        printf("%f\n",x[i]);

    printf("total number of iterations: %d\n", nit);
    cleanup();

    // Timing stuff
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("CPU Time used was: %f\n", cpu_time_used); // TODO remove after dealing with timing

    return EXIT_SUCCESS;
} // End of the main function
