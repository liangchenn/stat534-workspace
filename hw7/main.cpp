#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "matrices.h"

int main()
{
    const int n = 158;
    const int p = 51;
    const int num_draws = 10000;

    char datafilename[] = "erdata.txt";
    char sigma_output_file[] = "covariance.txt";
    char random_sample_output_file[] = "random-sample.txt";
    char random_sample_cov_output_file[] = "random-sample-covariance.txt";

    // problem 1 containers
    gsl_matrix *data = gsl_matrix_alloc(n, p);
    gsl_matrix *sigma = gsl_matrix_alloc(p, p);

    // problem 3 containers
    gsl_matrix *samples = gsl_matrix_alloc(num_draws, p);
    gsl_matrix *sample_cov = gsl_matrix_alloc(p, p);

    // --------------------------- problem 1: covariance matrix ---------------------------
    // load data
    load_data(datafilename, data);
    makeCovariance(sigma, data);
    save_matrix(sigma_output_file, sigma);

    // --------------------------- problem 2: multivariate normal RNG ---------------------------

    // 0. setups
    const gsl_rng_type *T = gsl_rng_mt19937;
    gsl_rng *rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, 5566); // seed

    // --------------------------- problem3: 10000 random samples covariance ---------------------------

    // 1. draw random samples and save
    randomMVN(rng, samples, sigma);
    save_matrix(random_sample_output_file, samples);

    makeCovariance(sample_cov, samples);
    save_matrix(random_sample_cov_output_file, sample_cov);

    // --------------------------- release mem ---------------------------
    gsl_rng_free(rng);
    gsl_matrix_free(data);
    gsl_matrix_free(sigma);
    gsl_matrix_free(samples);
    gsl_matrix_free(sample_cov);

    return 0;
}
