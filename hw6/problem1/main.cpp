#include <stdio.h>
#include <stdlib.h>

#include "matrices.h"

int main()
{
    // 1. load data
    int n = 10;
    int p = 10;
    char datafilename[] = "mybandedmatrix.txt";

    gsl_matrix *matrix = gsl_matrix_alloc(n, p);
    FILE *datafile = fopen(datafilename, "r");
    if (datafile == NULL)
    {
        printf("Error: could not open %s.\n", datafilename);
        return (1);
    }

    if (gsl_matrix_fscanf(datafile, matrix) != 0)
    {
        printf("Error: could not read matrix from %s.\n", datafilename);
        fclose(datafile);
        gsl_matrix_free(matrix);
        return (1);
    }
    fclose(datafile);

    // 2. calculate det of matrix
    double det = getDeterminant(matrix);
    printf("The determinant of %s is %.10lf.\n", datafilename, det);

    // 3. release memory
    gsl_matrix_free(matrix);

    return (0);
}
