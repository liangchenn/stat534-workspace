#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "regmodels.h"
#include "matrices.h"

int main()
{
    // 0. constants
    const int nMaxRegressions = 10;

    // 1. load data
    int i, j;

    int n = 158;                        // sample size
    int p = 51;                         // number of variables
    char datafilename[] = "erdata.txt"; // name of the data file
    char outputfilename[] = "top-10-regression-models.txt";
    char columnNameStorage[51][16];
    const char *columnNames[51];

    snprintf(columnNameStorage[0], sizeof(columnNameStorage[0]), "Y");
    columnNames[0] = columnNameStorage[0];
    for (i = 1; i < p; i++)
    {
        snprintf(columnNameStorage[i], sizeof(columnNameStorage[i]), "X%d", i);
        columnNames[i] = columnNameStorage[i];
    }

    // allocate the data matrix
    gsl_matrix *data = gsl_matrix_alloc(n, p);

    // read the data
    FILE *datafile = fopen(datafilename, "r");
    // readmatrix(datafilename, n, p, data);
    gsl_matrix_fscanf(datafile, data);
    fclose(datafile);

    // 2. Add regression model for sorted, non-dup. regressors

    // (0) max counter
    int counter = 0;

    // (1) create linked list
    RegressionPointer regressions = new Regression;

    regressions->lenA = 0;
    regressions->logmarglik = 0.0;
    regressions->A = NULL;
    regressions->Next = NULL;

    // (2) single predictor case
    int lenA = 1;
    int A[2];
    for (i = 1; i < p; i++) // first column is response
    {
        A[0] = i + 1;
        double logmarglik = marglik(data, lenA, A);
        AddRegression(nMaxRegressions, regressions, lenA, A, logmarglik, counter);
    }

    // (3) two-predictor case
    lenA = 2;

    for (i = 1; i < p; i++)
    {
        for (j = i + 1; j < p; j++)
        {
            // create non-duplicated index pair by only including {i, j > i}
            A[0] = i + 1;
            A[1] = j + 1;

            double logmarglik = marglik(data, lenA, A);
            AddRegression(nMaxRegressions, regressions, lenA, A, logmarglik, counter);
        }
    }

    // (4) output result
    SaveRegressions(outputfilename, regressions, columnNames);

    // (5) release mem
    DeleteAllRegressions(regressions, counter);
    gsl_matrix_free(data);

    return (0);
}
