#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

/* By Jason Yao */
typedef enum { false, true } bool; // Enables boolean types

/* Configuration settings */
bool IS_DEBUG_MODE = false;         /* Change to true to see intermediate values */
bool IS_SEQUENTIAL_MODE = false;    /* Change to true to switch to sequential version */
bool IS_TIMED_MODE = false;         /* Change to true to switch to show runtime */

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
        printf("%f ", curr[i]);
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
        printf("%f ", x[i]);
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
 * Parallel solution:
 * 1.) Based off of rank, if my_rank >= number_of_elements, then allocate extra spots based off of excessive threads
 * 2.) Build upper and lower ranges that each thread does. If my_rank >= number_of_elements, have the ranges be [0, 0].
 * 3.) If it is a valid working thread (defined as one with rank < num), have it calculate the x_i's in the given range
 * If it is a non-working thread, have it sit idle. (we set it idle due to no work able to be given).
 * 4.) Each one of the threads will then hit a syncronisation point in the form of MPI_Barrier.
 * 5.) After the barrier, each thread will then have their locally calculated x_i values copied and concatenated via
 * MPI_Allgather. For idle threads, they will copy the value "0", since apparently MPI_Allgather requires some data be
 * sent.
 * 6.) Each thread will then have the full new set of x_i, and can do their own computations of the errors values.
 * If all error values are less than the absolute value, then the thread finishes. If not, it repeats.
 * 7.) Only a single thread (arbitrarily 0 is chosen) will output the final information.
 */
int parallel2(int comm_size, int my_rank) {
    int count = 0;
    int local_n;

    // Calculates the size of the subarray that this process needs to deal with.
    if (comm_size > num)
    {
        // More processes than number of elements, give 1 element to each process with rank < num, and a white value
        local_n = 1;
    }
    else
    {
        // Less processes than number of elements, give ceil(number_of_elements/comm_size) to each process
        local_n = (num + comm_size - 1) / (comm_size); // ceil(number_of_elements/comm_size)
    }

    int specific_n = local_n;
    int invalidUseForN = local_n;

    // Indices to work between
    int lower_index_range = my_rank * local_n;
    int upper_index_range;

    // Checks if it's the last valid thread
    int numberOfValidProcesses = comm_size;
    if (comm_size > num)
        numberOfValidProcesses = num;

    // Identifies invalid threads
    int invalidThreadCount = 0;

    for (int i = 0; i < numberOfValidProcesses; ++i)
    {
        if (((i - 1) < num) && (i * invalidUseForN >= num))
        {
            ++invalidThreadCount;
            if (my_rank == i)
            {
                upper_index_range = 0;
                lower_index_range = 0;
                local_n = 1;
                specific_n = 1;
            }
        }
    }

    /**
     * Marks a process as invalid if they have nothing.
     * e.g. num = 5, p = 4: ceil(num/comm_sz) == 2
     * p_0: {A, B}
     * p_1: {C, D}
     * p_2: {E}
     * p_3: {} <---- We need to thus mark this thread as an invalid one.
     */
    if ((comm_size > num) && (invalidThreadCount == 1))
        numberOfValidProcesses = num - 1;
    else
        numberOfValidProcesses -= invalidThreadCount;

    if ((my_rank == numberOfValidProcesses - 1) && (num % numberOfValidProcesses != 0))
    {
        // Sets upper index range for last valid thread that isn't easily modulo'd
        upper_index_range = my_rank * local_n + (num % local_n);
        local_n = num % local_n;
    }
    else if ((my_rank == numberOfValidProcesses - 1) && (num % numberOfValidProcesses == 0)) {
        // Sets upper index range for last valid thread that is easily modulo'd
        upper_index_range = my_rank * local_n + local_n;
    }
    else if (my_rank >= numberOfValidProcesses)
    {
        // Sets upper index range for threads that don't do any work
        upper_index_range = 0;
    }
    else if (my_rank < numberOfValidProcesses)
    {
        // Sets upper index range for threads that do work
        upper_index_range = my_rank * local_n + local_n;
    }

    if (IS_DEBUG_MODE)
    {
        printf("my_rank: %d, local_n: %d, count: %d, specific_n: %d, "
                       "leftBound: %d, rightBound: %d, numberOfValidProcesses: %d\n", my_rank, local_n, count,
               specific_n, lower_index_range, upper_index_range, numberOfValidProcesses);
        printf("The number of valid threads: %d, Invalid threads: %d\n", numberOfValidProcesses, invalidThreadCount);
    }

    // Allocations
    float* local_new;
    int allocationAmount = num;
    if (comm_size > num)
        allocationAmount = comm_size;
    float* all_new = (float *) malloc(allocationAmount * sizeof(float));
    float* new_errors = (float *) malloc(allocationAmount * sizeof(float)); // For debugging only
    memset(new_errors, 0, sizeof(float) * allocationAmount);
    memset(all_new, 0, sizeof(float) * allocationAmount);

    // Reset of allocations
    if (my_rank < num)
        local_new = (float *) malloc(local_n * sizeof(float));
    else
        local_new = (float *) malloc(1 * sizeof(float));

    // Beginning of iteration
    for (; !SOLUTION_IS_SOLVED; ++count)
    {
        // Generates x_i if it is a process marked with work
        if (my_rank < num)
        {
            // Generates x_i if process is lower than comm_world_size
            int subarrayIndex = 0;
            for (int i = lower_index_range; i < upper_index_range; ++i)
            {
                // Finds x_i
                local_new[subarrayIndex] = b[i];

                for (int j = 0; j < num; ++j)
                    if (j != i)
                        local_new[subarrayIndex] -= (a[i][j] * x[j]);

                local_new[subarrayIndex] /= a[i][i];
                ++subarrayIndex;
            }
        } // End of x_i generation by working threads

        // Same block is reached by all processes
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgather(local_new, local_n, MPI_FLOAT, all_new, invalidUseForN, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // Each process checks each x_i for closeness validity, and checks whether it is complete for them
        bool isDone = true;

        for (int i = 0; i < num; ++i)
        {
            new_errors[i] = fabs((all_new[i] - x[i])/(all_new[i]));
            if (new_errors[i] > err)
                isDone = false;
            x[i] = all_new[i];
        }

        if ((IS_DEBUG_MODE) && (my_rank == 0))
            testGeneratedXAndErr(all_new, new_errors);

        if (isDone)
            SOLUTION_IS_SOLVED = true;
    } // End of iterating until done

    // Internal function cleanup
    free(new_errors);
    free(all_new);
    free(local_new);

    return count;
} // End of parallel2

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

    /* MPI stuff */
    int comm_size;      // Number of processes
    int my_rank;        // My process rank

    /* Solves the matrix */
    if (IS_SEQUENTIAL_MODE)
    {
        nit = sequential();
        /* Writing to the stdout */
        for(int i = 0; i < num; i++)
            printf("%f\n",x[i]);

        printf("total number of iterations: %d\n", nit);

        // Timing stuff
        if (IS_TIMED_MODE)
        {
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            printf("CPU Time used was: %f\n", cpu_time_used);
        }
    }
    else
    {
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        nit = parallel2(comm_size, my_rank);
        if (my_rank == 0)
        {
            /* Writing to the stdout */
            for(int i = 0; i < num; i++)
                printf("%f\n",x[i]);

            printf("total number of iterations: %d\n", nit);

            // Timing stuff
            if (IS_TIMED_MODE)
            {
                end = clock();
                cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
                printf("CPU Time used was: %f\n", cpu_time_used);
            }
        }
        MPI_Finalize();
    } // End of dealing with parallel

    cleanup();
    return EXIT_SUCCESS;
} // End of the main function
