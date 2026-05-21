// AI-generated response for Homework 5, Problem 1.
// Prompt: write double marglik(int n, int p, double** data, int lenA, int* A)
// for the marginal likelihood formula in hw5.pdf, using the provided matrix
// helper functions from Matrices/.

double marglik(int n, int p, double **data, int lenA, int *A)
{
    int i, j;
    double **XA = allocmatrix(n, lenA);
    double **XAT = NULL;
    double **MA = NULL;
    double **rhs = NULL;
    double **sol = NULL;

    (void)p;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < lenA; j++)
        {
            XA[i][j] = data[i][A[j] - 1];
        }
    }

    XAT = transposematrix(n, lenA, XA);
    MA = allocmatrix(lenA, lenA);
    matrixproduct(lenA, n, lenA, XAT, XA, MA);

    for (i = 0; i < lenA; i++)
    {
        MA[i][i] += 1.0;
    }

    double logdetMA = logdet(lenA, MA);
    double yty = 0.0;
    for (i = 0; i < n; i++)
    {
        yty += data[i][0] * data[i][0];
    }

    rhs = allocmatrix(lenA, 1);
    for (i = 0; i < lenA; i++)
    {
        for (j = 0; j < n; j++)
        {
            rhs[i][0] += XAT[i][j] * data[j][0];
        }
    }

    inverse(lenA, MA);
    sol = allocmatrix(lenA, 1);
    matrixproduct(lenA, lenA, 1, MA, rhs, sol);

    double quad = 0.0;
    for (i = 0; i < lenA; i++)
    {
        quad += rhs[i][0] * sol[i][0];
    }

    double alpha0 = (lenA + 2.0) / 2.0;
    double alphaN = (n + lenA + 2.0) / 2.0;
    double scale = 1.0 + yty - quad;
    double out = lgamma(alphaN) - lgamma(alpha0) - 0.5 * logdetMA
                 - alphaN * log(scale);

    freematrix(n, XA);
    freematrix(lenA, XAT);
    freematrix(lenA, MA);
    freematrix(lenA, rhs);
    freematrix(lenA, sol);

    return out;
}
