#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

void ai_makeCovariance(gsl_matrix *covX, gsl_matrix *X)
{
    size_t n = X->size1;
    size_t p = X->size2;

    if (covX->size1 != p || covX->size2 != p)
    {
        fprintf(stderr, "Error: covariance matrix has wrong dimension.\n");
        exit(1);
    }

    double *means = new double[p];

    for (size_t j = 0; j < p; j++)
    {
        double sum = 0.0;
        for (size_t i = 0; i < n; i++)
        {
            sum += gsl_matrix_get(X, i, j);
        }
        means[j] = sum / n;
    }

    for (size_t j = 0; j < p; j++)
    {
        for (size_t k = 0; k < p; k++)
        {
            double cov = 0.0;
            for (size_t i = 0; i < n; i++)
            {
                cov += (gsl_matrix_get(X, i, j) - means[j]) *
                       (gsl_matrix_get(X, i, k) - means[k]);
            }
            gsl_matrix_set(covX, j, k, cov / (n - 1));
        }
    }

    delete[] means;
}

gsl_matrix *ai_makeCholesky(gsl_matrix *K)
{
    if (K->size1 != K->size2)
    {
        fprintf(stderr, "Error: Cholesky decomposition needs a square matrix.\n");
        exit(1);
    }

    gsl_matrix *L = gsl_matrix_alloc(K->size1, K->size2);
    gsl_matrix_memcpy(L, K);

    int status = gsl_linalg_cholesky_decomp(L);
    if (status != 0)
    {
        fprintf(stderr, "Error: Cholesky decomposition failed.\n");
        gsl_matrix_free(L);
        exit(1);
    }

    for (size_t i = 0; i < L->size1; i++)
    {
        for (size_t j = i + 1; j < L->size2; j++)
        {
            gsl_matrix_set(L, i, j, 0.0);
        }
    }

    return L;
}

void ai_randomMVN(gsl_rng *mystream, gsl_matrix *samples, gsl_matrix *sigma)
{
    size_t num_samples = samples->size1;
    size_t p = sigma->size1;

    if (samples->size2 != p || sigma->size2 != p)
    {
        fprintf(stderr, "Error: sample and covariance dimensions do not match.\n");
        exit(1);
    }

    gsl_matrix *L = ai_makeCholesky(sigma);
    gsl_vector *z = gsl_vector_alloc(p);
    gsl_vector *x = gsl_vector_alloc(p);

    for (size_t i = 0; i < num_samples; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            gsl_vector_set(z, j, gsl_ran_ugaussian(mystream));
        }

        gsl_blas_dgemv(CblasNoTrans, 1.0, L, z, 0.0, x);

        for (size_t j = 0; j < p; j++)
        {
            gsl_matrix_set(samples, i, j, gsl_vector_get(x, j));
        }
    }

    gsl_vector_free(z);
    gsl_vector_free(x);
    gsl_matrix_free(L);
}

void ai_makeRandomSampleCovariance(gsl_rng *mystream,
                                   gsl_matrix *sample_cov,
                                   gsl_matrix *sigma,
                                   int num_draws)
{
    gsl_matrix *samples = gsl_matrix_alloc(num_draws, sigma->size1);

    ai_randomMVN(mystream, samples, sigma);
    ai_makeCovariance(sample_cov, samples);

    gsl_matrix_free(samples);
}
