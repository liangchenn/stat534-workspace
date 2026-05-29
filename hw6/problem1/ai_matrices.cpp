#include <stdio.h>
#include <stdlib.h>

#include "ai_matrices.h"

double ai_getDeterminant(gsl_matrix *matrix)
{
    int n = matrix->size1;

    if (matrix->size1 != matrix->size2)
    {
        printf("Error: determinant needs a square matrix.\n");
        exit(1);
    }

    if (n == 1)
    {
        return gsl_matrix_get(matrix, 0, 0);
    }

    if (n == 2)
    {
        return gsl_matrix_get(matrix, 0, 0) * gsl_matrix_get(matrix, 1, 1) -
               gsl_matrix_get(matrix, 0, 1) * gsl_matrix_get(matrix, 1, 0);
    }

    double det = 0.0;

    for (int col = 0; col < n; col++)
    {
        double value = gsl_matrix_get(matrix, 0, col);
        if (value == 0.0)
        {
            continue;
        }

        gsl_matrix *minor = ai_getSubmatrix(matrix, 0, col);
        double sign = (col % 2 == 0) ? 1.0 : -1.0;
        det += sign * value * ai_getDeterminant(minor);
        gsl_matrix_free(minor);
    }

    return det;
}

gsl_matrix *ai_getSubmatrix(gsl_matrix *matrix, int removeRow, int removeCol)
{
    int n = matrix->size1;
    gsl_matrix *minor = gsl_matrix_alloc(n - 1, n - 1);
    int outRow = 0;

    for (int i = 0; i < n; i++)
    {
        if (i == removeRow)
        {
            continue;
        }

        int outCol = 0;
        for (int j = 0; j < n; j++)
        {
            if (j == removeCol)
            {
                continue;
            }

            gsl_matrix_set(minor, outRow, outCol, gsl_matrix_get(matrix, i, j));
            outCol++;
        }

        outRow++;
    }

    return minor;
}
