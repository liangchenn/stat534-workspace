// AI-generated response for Homework 5, Problem 2.
// Prompt: write double marglik(gsl_matrix* data, int lenA, int* A)
// for the marginal likelihood formula in hw5.pdf, using GSL matrix routines
// and the helper functions from MatricesGSL/.

double marglik(gsl_matrix *data, int lenA, int *A)
{
    int i, j;
    int n = data->size1;
    gsl_matrix *XA = gsl_matrix_alloc(n, lenA);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < lenA; j++)
        {
            gsl_matrix_set(XA, i, j, gsl_matrix_get(data, i, A[j] - 1));
        }
    }

    gsl_matrix *XAT = transposematrix(XA);
    gsl_matrix *MA = gsl_matrix_alloc(lenA, lenA);
    matrixproduct(XAT, XA, MA);

    for (i = 0; i < lenA; i++)
    {
        gsl_matrix_set(MA, i, i, gsl_matrix_get(MA, i, i) + 1.0);
    }

    double logdetMA = logdet(MA);
    double yty = 0.0;
    for (i = 0; i < n; i++)
    {
        double y = gsl_matrix_get(data, i, 0);
        yty += y * y;
    }

    gsl_matrix *rhs = gsl_matrix_alloc(lenA, 1);
    for (i = 0; i < lenA; i++)
    {
        double value = 0.0;
        for (j = 0; j < n; j++)
        {
            value += gsl_matrix_get(XAT, i, j) * gsl_matrix_get(data, j, 0);
        }
        gsl_matrix_set(rhs, i, 0, value);
    }

    gsl_matrix *MAinv = inverse(MA);
    gsl_matrix *sol = gsl_matrix_alloc(lenA, 1);
    matrixproduct(MAinv, rhs, sol);

    double quad = 0.0;
    for (i = 0; i < lenA; i++)
    {
        quad += gsl_matrix_get(rhs, i, 0) * gsl_matrix_get(sol, i, 0);
    }

    double alpha0 = (lenA + 2.0) / 2.0;
    double alphaN = (n + lenA + 2.0) / 2.0;
    double scale = 1.0 + yty - quad;
    double out = lgamma(alphaN) - lgamma(alpha0) - 0.5 * logdetMA
                 - alphaN * log(scale);

    gsl_matrix_free(XA);
    gsl_matrix_free(XAT);
    gsl_matrix_free(MA);
    gsl_matrix_free(rhs);
    gsl_matrix_free(MAinv);
    gsl_matrix_free(sol);

    return out;
}
