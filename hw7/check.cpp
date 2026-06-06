#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_matrix.h>

int load_matrix(char *filename, gsl_matrix *matrix)
{
    FILE *infile = fopen(filename, "r");
    if (infile == NULL)
    {
        fprintf(stderr, "Error: could not open %s.\n", filename);
        return 1;
    }

    if (gsl_matrix_fscanf(infile, matrix) != 0)
    {
        fprintf(stderr, "Error: could not read matrix from %s.\n", filename);
        fclose(infile);
        return 1;
    }

    fclose(infile);
    return 0;
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

int main()
{
    const int p = 51;

    char sigma_file[] = "covariance.txt";
    char sample_cov_file[] = "random-sample-covariance.txt";
    char abs_diff_file[] = "absolute-difference.txt";

    gsl_matrix *sigma = gsl_matrix_alloc(p, p);
    gsl_matrix *sample_cov = gsl_matrix_alloc(p, p);
    gsl_matrix *abs_diff = gsl_matrix_alloc(p, p);

    if (sigma == NULL || sample_cov == NULL || abs_diff == NULL)
    {
        fprintf(stderr, "Error: failed to allocate matrices.\n");
        gsl_matrix_free(sigma);
        gsl_matrix_free(sample_cov);
        gsl_matrix_free(abs_diff);
        return 1;
    }

    if (load_matrix(sigma_file, sigma) != 0 ||
        load_matrix(sample_cov_file, sample_cov) != 0)
    {
        gsl_matrix_free(sigma);
        gsl_matrix_free(sample_cov);
        gsl_matrix_free(abs_diff);
        return 1;
    }

    double sum_abs_diff = 0.0;
    double max_abs_diff = 0.0;
    size_t max_i = 0;
    size_t max_j = 0;

    for (size_t i = 0; i < p; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            double diff = fabs(gsl_matrix_get(sigma, i, j) -
                               gsl_matrix_get(sample_cov, i, j));

            gsl_matrix_set(abs_diff, i, j, diff);
            sum_abs_diff += diff;

            if (diff > max_abs_diff)
            {
                max_abs_diff = diff;
                max_i = i;
                max_j = j;
            }
        }
    }

    if (save_matrix(abs_diff_file, abs_diff) != 0)
    {
        gsl_matrix_free(sigma);
        gsl_matrix_free(sample_cov);
        gsl_matrix_free(abs_diff);
        return 1;
    }

    printf("Mean absolute difference: %.6f\n", sum_abs_diff / (p * p));
    printf("Max absolute difference:  %.6f\n", max_abs_diff);
    printf("Max location: row %zu, col %zu\n", max_i + 1, max_j + 1);
    printf("Saved absolute difference matrix to %s\n", abs_diff_file);

    gsl_matrix_free(sigma);
    gsl_matrix_free(sample_cov);
    gsl_matrix_free(abs_diff);

    return 0;
}
