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

    if ((comm_size > num) && (my_rank < num))
        local_n = 1;
    else if ((comm_size > num) && (my_rank >= num))
        local_n = 1;
    else
        local_n = (num + comm_size - 1) / (comm_size); // sub_array_size_ceiling, number of elements in the subarray

    int specific_n = local_n;

    // Indices to work between
    int lower_index_range = my_rank * local_n;
    int upper_index_range;

    // Checks if it's the last valid thread
    int numberOfValidProcesses = comm_size;
    if (comm_size > num)
        numberOfValidProcesses = num;

    if ((my_rank == numberOfValidProcesses - 1) && (num % numberOfValidProcesses != 0)) {
//        printf("GOT SET HERE\n"); //TODO remove after
        // Sets upper index range for last valid thread that isn't easily modulo'd
        upper_index_range = my_rank * local_n + (num % local_n);
        local_n = num % local_n;
    }
    else if ((my_rank == numberOfValidProcesses - 1) && (num % numberOfValidProcesses == 0)) {
        // Sets upper index range for last valid thread that is easily modulo'd
        upper_index_range = my_rank * local_n + local_n;
    }
    else if (my_rank >= num)
    {
        // Sets upper index range for threads that don't do any work
        upper_index_range = 0;
    }
    else if (my_rank < numberOfValidProcesses - 1)
    {
        // Sets upper index range for threads that do work
        upper_index_range = my_rank * local_n + local_n;
    }

    //TODO remove after
//    if (my_rank >= num)
//    {
//        printf("my_rank: %d\n", my_rank);
//        printf("local_n: %d, count: %d, specific_n: %d, leftBound: %d, rightBound: %d, numberOfValidProcesses: %d\n",
//               local_n, count, specific_n, lower_index_range, upper_index_range, numberOfValidProcesses);
//    }

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

//        printf("Check 1: my_rank = %d\n", my_rank); //TODO remove after
        // Same block is reached by all processes
        MPI_Barrier(MPI_COMM_WORLD);
//        printf("Check 2: my_rank = %d\n", my_rank); // TODO remove after
        MPI_Allgather(local_new, local_n, MPI_FLOAT, all_new, specific_n, MPI_FLOAT, MPI_COMM_WORLD);
//        printf("Check 3: my_rank = %d\n", my_rank); //TODO remove after
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

    free(new_errors);
    free(all_new);
    free(local_new);

    return count;
} // End of parallel2

int parallel(int comm_size, int my_rank)
{

    int count = 0;
    int local_n;

    if (comm_size > num)
        local_n = 1;
    else
        local_n = (num + comm_size - 1)/(comm_size); // sub_array_size_ceiling, number of elements in the subarray

    int specific_n = local_n;

    // Indices to work between
    int lower_index_range = my_rank * local_n;
    int upper_index_range;

    // Checks if it's the last valid thread
    int numberOfValidProcesses = comm_size;
    if (comm_size > num)
        numberOfValidProcesses = num;

    if ((my_rank == numberOfValidProcesses - 1) && (num % numberOfValidProcesses != 0))
    {
        // Sets upper index range for last valid thread that isn't easily modulo'd
        upper_index_range = my_rank * local_n + (num % local_n);
        local_n = num % local_n;
    }
    else
        upper_index_range = my_rank * local_n + local_n;

    float* local_new = (float *) malloc(local_n * sizeof(float));
    float* all_new = (float *) malloc(num * sizeof(float));
    float* new_errors = (float *) malloc(num * sizeof(float)); // For debugging only
    memset(new_errors, 0, sizeof(float) * num);
    memset(all_new, 0, sizeof(float) * num);

    // Beginning of iteration
    for (; !SOLUTION_IS_SOLVED; ++count)
    {
        // Generates all new x_i for the given subarray
        int subarrayIndex = 0;
        for (int i = lower_index_range; i < upper_index_range; ++i)
        {
            // Finds x_i
            local_new[subarrayIndex] = b[i];
            printf("------------------------ local_new[subarrayIndex] = %f\n", local_new[subarrayIndex]);

            for (int j = 0; j < num; ++j)
                if (j != i)
                    local_new[subarrayIndex] -= (a[i][j] * x[j]);

            local_new[subarrayIndex] /= a[i][i];
            ++subarrayIndex;
        }

        // Gathers each processes's local_new[], concatenates, and stores in each process's all_new[]
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgather(local_new, local_n, MPI_FLOAT, all_new, specific_n, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

//        //TODO remove after testing
//        if (my_rank == comm_size - 1)
//        {
//            printf("\n-------------------BIG TEST------------------------------\n");
//            for (int i = 0; i < num; ++i)
//                printf("%f ", all_new[i]);
//            printf("\n");
//            printf("\n-------------------BIG TEST------------------------------\n");
//        }

        // Checks each x_i for closeness validity
        bool isDone = true;
        for (int i = 0; i < num; ++i)
        {
            new_errors[i] = fabs((all_new[i] - x[i])/(all_new[i]));
            if (new_errors[i] > err)
                isDone = false;
            x[i] = all_new[i];
        }

        //TODO remove after testing: the thing is the indices are weird after completion
//        if (my_rank == comm_size - 1)
//        {
//            printf("\n\n\n\n");
//            for (int i = 0; i < num; ++i)
//                printf("x[%d] is: %f\n", i, x[i]);
//            printf("\n\n");
//            for (int i = 0; i < num; ++i)
//                printf("local_new[%d] is: %f\n", i, local_new[i]);
//            printf("\n\n\n\n");
//        }

        if ((IS_DEBUG_MODE) && (my_rank == 0))
            testGeneratedXAndErr(all_new, new_errors);

        if (isDone)
            SOLUTION_IS_SOLVED = true;
    } // End of iterating until completed

    free(new_errors);
    free(all_new);
    free(local_new);

    return count;
} // End of parallel 2



//int parallel(int comm_size, int my_rank)
//{
//    int count = 0;
//    int local_n = (num + comm_size - 1)/(comm_size); // sub_array_size_ceiling, number of elements in the subarray
//    int specific_n = local_n;
//
//    // Indices to work between
//    int lower_index_range = my_rank * local_n;
//    int upper_index_range;
//
//    // Checks if it's the last one
//    if ((my_rank == comm_size - 1) && (num % comm_size != 0))
//    {
//        upper_index_range = my_rank * local_n + (num % local_n);
//        local_n = num % local_n;
//    }
//    else
//        upper_index_range = my_rank * local_n + local_n;
//
//    float* local_new = (float *) malloc(local_n * sizeof(float));
//    float* all_new = (float *) malloc(num * sizeof(float));
//    float* new_errors = (float *) malloc(num * sizeof(float)); // For debugging only
//    memset(new_errors, 0, sizeof(float) * num);
//    memset(all_new, 0, sizeof(float) * num);
//
////    printf("comm_size is: %d, number of elements is: %d\n", comm_size, num); //TODO remove after
//
//    // Beginning of iteration
//    for (; !SOLUTION_IS_SOLVED; ++count)
//    {
////        // TODO remove after
////        if (count == 2)
////            break;
//
//        // Generates all new x_i for the given subarray
//        int subarrayIndex = 0;
//        for (int i = lower_index_range; i < upper_index_range; ++i)
//        {
//            // Finds x_i
//            local_new[subarrayIndex] = b[i];
//
//            //TODO remove after
//            // Checks b[i] and local_new[subArrayIndex]
////            if ((my_rank == comm_size - 1) && (i == 8 || i == 9))
////                printf("####################### b[%d] = %f, local_new[subarray] = %f\n", i, b[i], local_new[subarrayIndex]);
//
//            for (int j = 0; j < num; ++j)
//                if (j != i)
//                {
//                    local_new[subarrayIndex] -= (a[i][j] * x[j]);
//
////                    //TODO remove after
////                    if ((my_rank == comm_size - 1) && (i == 8 || i == 9))
////                    {
////                        printf("i: %d, j: %d, Local_new[%d] = %f\n", i, j, i, local_new[subarrayIndex]);
////                        printf("a[%d][%d] is: %f\n", i, j, a[i][j]);
////                        printf("x[%d] is: %f\n", j, x[j]);
////                    }
//                }
//
//            //printf("$$$$$$$$$$$$$$$$ local_new[%d] is: %f\n", i, local_new[subarrayIndex]);
//
////            local_new[subarrayIndex] /= a[i][i];
////            ++subarrayIndex;
//
//            //TODO remove after testing
////            if ((my_rank == comm_size - 1) && (i == 8 || i == 9))
////            {
////                printf("\n-------------------------------------------------\n");
////                printf("Subarray index = %d\n", subarrayIndex);
////                printf("my_rank is: %d, local_n is: %d\n", my_rank, local_n);
////                printf("index %d is complete, value calculated was: %f", i, local_new[subarrayIndex - 1]);
////                printf("\n-------------------------------------------------\n");
////            }
//        }
//
//        // Gathers each processes's local_new[], concatenates, and stores in each process's all_new[]
////        MPI_Barrier(MPI_COMM_WORLD);
////        MPI_Allgather(local_new, local_n, MPI_FLOAT, all_new, specific_n, MPI_FLOAT, MPI_COMM_WORLD);
//        //MPI_Barrier(MPI_COMM_WORLD);
//
////        //TODO remove after testing
////        if (my_rank == comm_size - 1)
////        {
////            printf("\n-------------------BIG TEST------------------------------\n");
////            for (int i = 0; i < num; ++i)
////                printf("%f ", all_new[i]);
////            printf("\n");
////            printf("\n-------------------BIG TEST------------------------------\n");
////        }
//
//        // Checks each x_i for closeness validity
////        bool isDone = true;
////        for (int i = 0; i < num; ++i)
////        {
////            new_errors[i] = fabs((all_new[i] - x[i])/(all_new[i]));
////            if (new_errors[i] > err)
////                isDone = false;
////            x[i] = all_new[i];
////        }
//
//        //TODO remove after testing: the thing is the indices are weird after completion
////        if (my_rank == comm_size - 1)
////        {
////            printf("\n\n\n\n");
////            for (int i = 0; i < num; ++i)
////                printf("x[%d] is: %f\n", i, x[i]);
////            printf("\n\n");
////            for (int i = 0; i < num; ++i)
////                printf("local_new[%d] is: %f\n", i, local_new[i]);
////            printf("\n\n\n\n");
////        }
//
//        if ((IS_DEBUG_MODE) && (my_rank == 0))
//            testGeneratedXAndErr(all_new, new_errors);
//
////        if (isDone)
//        SOLUTION_IS_SOLVED = true;
//    } // End of iterating until completed
//
//    free(new_errors);
//    free(all_new);
//    free(local_new);
//
//    return count;
//} // End of parallel 2


/**
 * Parallel version of the code to solve n equations to within a bounds
 */
int oldParallel(int comm_size, int my_rank)
{
//    printf("#################\n");
//    // Sets the current x values to the initial x values
//    for (int i = 0; i < num; ++i)
//        curr[i] = x[i];
//    testCurrent();
//    printf("#################\n");

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

    int count = 0;
    int local_n = (num + comm_size - 1)/(comm_size); // sub_array_size_ceiling, number of elements in the subarray
    float* local_curr = (float *) malloc(local_n * sizeof(float)); // recieving buffer for elements being scattered in blocks
    float* local_x = (float *) malloc(local_n * sizeof(float)); // The newly calculated values locally
    float* new_x = (float *) malloc(num * sizeof(float));       // Contains all new values at p0
    float* new_err = (float *) malloc(num * sizeof(float));     // Contains all new errors

    printf("Local_n is: %d\n", local_n);
    for (; !SOLUTION_IS_SOLVED; ++count)
    {
        // Scatters the current x_i values to each process
        MPI_Scatter(curr, local_n, MPI_FLOAT, local_curr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Generates all new x_i
        for (int i = 0; i < local_n; ++i)
        {
            // Finds x_i
            int actualIndex = (local_n * my_rank) + i;
            if (actualIndex > num)
                break;

            local_x[i] = b[actualIndex];

            for (int j = 0; j < num; ++j)
                if (j != i)
                    local_x[i] -= (a[i][j] * curr[j]);
            local_x[i] /= a[i][i];
        } // At the end of this, local_x[] holds the newly calculated x_i subarray

        // Gathers each processes's local_x[], concatenates, and stores in each process's curr[]
        MPI_Allgather(local_x, local_n, MPI_FLOAT, new_x, local_n, MPI_FLOAT, MPI_COMM_WORLD);

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
    } // End of iterating until done

    /* Working */
//    if (my_rank == 0)
//    {
//        // Process is root
//        for (; !SOLUTION_IS_SOLVED; ++count)
//        {
//            // Scatters the current x_i values to each process
//            MPI_Scatter(curr, local_n, MPI_FLOAT, local_curr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
//
//            // Generates all new x_i
//            for (int i = 0; i < local_n; ++i)
//            {
//                // Finds x_i
//                local_x[i] = b[(local_n * my_rank) + i];
//
//                for (int j = 0; j < num; ++j)
//                    if (j != i)
//                        local_x[i] -= (a[i][j] * curr[j]);
//                local_x[i] /= a[i][i];
//            } // At the end of this, local_x[] holds the newly calculated x_i subarray
//
//            // Gathers each processes's local_x[], concatenates, and stores in each process's curr[]
//            MPI_Allgather(local_x, local_n, MPI_FLOAT, new_x, local_n, MPI_FLOAT, MPI_COMM_WORLD);
//
//            // Checks each x_i for closeness validity
//            bool isDone = true;
//            for (int i = 0; i < num; ++i)
//            {
//                new_err[i] = fabs((new_x[i] - curr[i])/(new_x[i]));
//                if (new_err[i] > err)
//                    isDone = false;
//
//                curr[i] = new_x[i];
//                x[i] = curr[i];
//            }
//
//            if (IS_DEBUG_MODE)
//                testGeneratedXAndErr(new_x, new_err);
//
//            if (isDone)
//                SOLUTION_IS_SOLVED = true;
//        } // End of iterating until done
//    } // End of dealing with root rank
//    else
//    {
//        for (; !SOLUTION_IS_SOLVED; ++count)
//        {
//            // Process is not root
//            MPI_Scatter(a, local_n, MPI_FLOAT, local_curr, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
//
//            // Generates all new x_i
//            for (int i = 0; i < local_n; ++i)
//            {
//                // Finds x_i
//                local_x[i] = b[(local_n * my_rank) + i];
//
//                for (int j = 0; j < num; ++j)
//                    if (j != i)
//                        local_x[i] -= (a[i][j] * curr[j]);
//                local_x[i] /= a[i][i];
//            } // At the end of this, local_x[] holds the newly calculated x_i subarray
//
//            // Gathers each processes's local_x[], concatenates, and stores in each process's curr[]
//            MPI_Allgather(local_x, local_n, MPI_FLOAT, new_x, local_n, MPI_FLOAT, MPI_COMM_WORLD);
//
//            // Checks each x_i for closeness validity
//            bool isDone = true;
//            for (int i = 0; i < num; ++i)
//            {
//                new_err[i] = fabs((new_x[i] - curr[i])/(new_x[i]));
//                if (new_err[i] > err)
//                    isDone = false;
//
//                curr[i] = new_x[i];
//                x[i] = curr[i];
//            }
//
//            if (IS_DEBUG_MODE)
//                testGeneratedXAndErr(new_x, new_err);
//
//            if (isDone)
//                SOLUTION_IS_SOLVED = true;
//        } // End of iterating until done
//    } // End of dealing with non-root rank

    // Frees everything
    free(local_curr);
    free(new_x);
    free(new_err);

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
    }
    else
    {
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


        /**
         * TODO NOTE: cannot do this, since everything needs to get to the same sync point, so a deadlock appears.
         * wrap generation in if statement, and have them sent to each
         */
        // Sanity check to kill off excess processes
//        printf("my_rank = %d\n", my_rank);
//        printf("num is %d\n", num);
//        if (my_rank >= num)
//        {
//            // Exits
//            printf("Process %d has exited early correctly\n", my_rank);
//            MPI_Finalize();
//            cleanup();
//            return EXIT_SUCCESS;
//        }

        printf("----------------------START: my_rank = %d\n", my_rank);

        nit = parallel2(comm_size, my_rank);

        printf("----------------------END: my_rank = %d\n", my_rank);

        if (my_rank == 0)
        {
            /* Writing to the stdout */
            for(int i = 0; i < num; i++)
                printf("%f\n",x[i]);

            printf("total number of iterations: %d\n", nit);

            // Timing stuff
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            printf("CPU Time used was: %f\n", cpu_time_used); // TODO remove after dealing with timing
        }
        MPI_Finalize();
    } // End of dealing with parallel

    cleanup();
    return EXIT_SUCCESS;
} // End of the main function
