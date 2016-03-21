#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* By Jason Yao */
typedef enum { false, true } bool; // Enables boolean types

/* Configuration settings */
bool IS_DEBUG_MODE = false;         /* Change to true to see intermediate values */
bool IS_SEQUENTIAL_MODE = true;     /* Change to true to switch to sequential version */

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
} // End of the generate next round function

/**
 * Parallel version of the code to solve n equations to within a bounds
 */
int parallel()
{
    //TODO do parallel version
    return 0;
}

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
} // End of the solve function

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
}


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


int main(int argc, char *argv[])
{

    int i;
    int nit = 0; /* number of iterations */

    if( argc != 2)
    {
        printf("Usage: gsref filename\n");
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
    for( i = 0; i < num; i++)
        printf("%f\n",x[i]);

    printf("total number of iterations: %d\n", nit);
    exit(0);
} // End of the main function
