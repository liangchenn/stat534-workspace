#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrices.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

int load_data(char *filename, gsl_matrix *data)
{
    FILE *datafile = fopen(filename, "r");
    if (datafile == NULL)
    {
        fprintf(stderr, "Error: could not open %s.\n", filename);
        return 1;
    }

    if (gsl_matrix_fscanf(datafile, data) != 0)
    {
        fprintf(stderr, "Error: could not read data matrix from %s.\n", filename);
        fclose(datafile);
        return 1;
    }

    fclose(datafile);
    return 0;
}

void makeCovariance(gsl_matrix *covX, gsl_matrix *X)
{
    size_t n = X->size1;
    size_t p = X->size2;

    // since covariance matrix is symmetric, (j, k) = (k, j), for j != k
    // diagonal case
    for (size_t j = 0; j < p; j++)
    {
        double cov_jj = gsl_stats_variance(gsl_matrix_const_ptr(X, 0, j), X->tda, n);
        gsl_matrix_set(covX, j, j, cov_jj);
    }
    // off-diagonal
    for (size_t j = 0; j < p; j++)
    {
        for (size_t k = j + 1; k < p; k++)
        {
            double cov_jk = gsl_stats_covariance(
                gsl_matrix_const_ptr(X, 0, j), X->tda,
                gsl_matrix_const_ptr(X, 0, k), X->tda,
                n);
            gsl_matrix_set(covX, j, k, cov_jk);
            gsl_matrix_set(covX, k, j, cov_jk);
        }
    }
}

gsl_matrix *makeCholesky(gsl_matrix *K)
{
    gsl_matrix *chol = gsl_matrix_alloc(K->size1, K->size2);
    gsl_matrix_memcpy(chol, K);

    gsl_linalg_cholesky_decomp(chol); // this modifies the input in-place
    // based on the documentation, it will only modify the lower triangular part
    // need to ignore the upper tri part when doing calculations afterwards.
    for (size_t i = 0; i < chol->size1; i++)
    {
        for (size_t j = i + 1; j < chol->size2; j++)
        {
            gsl_matrix_set(chol, i, j, 0.0);
        }
    }

    return chol;
}

void draw_multivariate_standard_normal(gsl_rng *rng, gsl_matrix *z)
{
    const size_t p = z->size1;

    for (size_t j = 0; j < p; j++)
    {
        const double zj = gsl_ran_ugaussian(rng);
        gsl_matrix_set(z, j, 0, zj);
    }
}

void randomMVN(gsl_rng *mystream, gsl_matrix *samples, gsl_matrix *sigma)
{
    size_t p = sigma->size1;
    gsl_matrix *chol = makeCholesky(sigma);
    gsl_matrix *Z = gsl_matrix_alloc(p, 1);
    gsl_matrix *X = gsl_matrix_alloc(p, 1);

    for (size_t i = 0; i < samples->size1; i++)
    {
        draw_multivariate_standard_normal(mystream, Z);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, chol, Z, 0.0, X);

        for (size_t j = 0; j < samples->size2; j++)
        {
            gsl_matrix_set(samples, i, j, gsl_matrix_get(X, j, 0));
        }
    }

    gsl_matrix_free(chol);
    gsl_matrix_free(Z);
    gsl_matrix_free(X);
}

int save_matrix(char *filename, gsl_matrix *matrix)
{
    FILE *outfile = fopen(filename, "w");
    if (outfile == NULL)
    {
        fprintf(stderr, "Error: could not open %s for writing.\n", filename);
        return 1;
    }

    for (size_t i = 0; i < matrix->size1; i++)
    {
        for (size_t j = 0; j < matrix->size2; j++)
        {
            fprintf(outfile, "%.6f%s", gsl_matrix_get(matrix, i, j),
                    (j + 1 == matrix->size2) ? "" : " ");
        }
        fprintf(outfile, "\n");
    }

    fclose(outfile);
    return 0;
}

// void print_comparison(gsl_matrix *sigma, gsl_matrix *sim_cov)
// {
//     double max_abs_diff = 0.0;
//     double sum_abs_diff = 0.0;
//     size_t count = sigma->size1 * sigma->size2;

//     for (size_t i = 0; i < sigma->size1; i++)
//     {
//         for (size_t j = 0; j < sigma->size2; j++)
//         {
//             double diff = fabs(gsl_matrix_get(sigma, i, j) -
//                                gsl_matrix_get(sim_cov, i, j));
//             sum_abs_diff += diff;
//             if (diff > max_abs_diff)
//             {
//                 max_abs_diff = diff;
//             }
//         }
//     }

//     printf("\nComparison between Sigma and simulated covariance matrix:\n");
//     printf("Mean absolute difference: %.6f\n", sum_abs_diff / count);
//     printf("Max absolute difference:  %.6f\n", max_abs_diff);
// }
