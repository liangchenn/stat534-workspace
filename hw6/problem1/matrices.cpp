#include <stdio.h>
#include <stdlib.h>

#include "matrices.h"

double getDeterminant(gsl_matrix *matrix)
{
    double det = 0.0;
    int n_rows = matrix->size1;
    int n_cols = matrix->size2;

    // basic check for squared matrix
    if (n_rows != n_cols)
    {
        printf("Error: matrix must be square, currently got (%d, %d).\n", n_rows, n_cols);
        exit(1);
    }

    // base case
    if (n_rows == 1)
    {
        return gsl_matrix_get(matrix, 0, 0);
    }

    // apply Laplace decomposition along i = 0 axis
    const int i = 0;

    for (int j = 0; j < n_cols; j++)
    {
        double aij = gsl_matrix_get(matrix, i, j);
        if (aij == 0)
        {
            // skip aij=0 case, no impact on determinant
            continue;
        }

        gsl_matrix *submatrix = getSubmatrix(matrix, i, j);
        double subdet = getDeterminant(submatrix);
        det += ((j % 2 == 0 ? 1.0 : -1.0) * aij * subdet);
        gsl_matrix_free(submatrix);
    }

    return det;
}

gsl_matrix *getSubmatrix(gsl_matrix *matrix, int row_remove_index, int col_remove_index)
{
    int n_rows = matrix->size1;
    int n_cols = matrix->size2;

    // matrix format check
    if (n_rows != n_cols)
    {
        printf("Error: matrix must be square, currently got (%d, %d).\n", n_rows, n_cols);
        exit(1);
    }

    // submatrix construction
    gsl_matrix *submatrix = gsl_matrix_alloc(n_rows - 1, n_cols - 1);

    int sub_i = 0;
    for (int i = 0; i < n_rows; i++)
    {
        if (i == row_remove_index)
        {
            // skip i-th row
            continue;
        }

        int sub_j = 0;
        for (int j = 0; j < n_cols; j++)
        {
            if (j == col_remove_index)
            {
                // skip j-th col
                continue;
            }

            // insert values for submatrix
            double value = gsl_matrix_get(matrix, i, j);
            gsl_matrix_set(submatrix, sub_i, sub_j, value);

            sub_j++;
        }

        sub_i++;
    }

    return submatrix;
}
