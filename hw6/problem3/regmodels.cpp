#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "regmodels.h"

void AddRegression(
    int nMaxRegs,
    RegressionPointer regression,
    int lenA,
    int *A,
    double logmarglik,
    int &counter)
{
    // 1. early stop case
    if (counter >= nMaxRegs &&
        regression->Next != NULL &&
        logmarglik <= regression->Next->logmarglik)
    {
        return;
    }

    // 2. insert position
    int i;
    RegressionPointer p = regression;
    RegressionPointer pnext = p->Next;

    while (pnext != NULL)
    {
        // TODO: check if regressor set A is finished, if true, return
        // NOTE: this is currently handled by creating sorted, non-duplicated index set in main.cpp

        // case: if new logmarglik is larger, go to next
        if (logmarglik > pnext->logmarglik)
        {
            p = pnext;
            pnext = p->Next;
        }
        else
        {
            // otherwise, stop and insert node
            break;
        }
    }

    // 3. create new regression node
    RegressionPointer pnew = new Regression;
    pnew->lenA = lenA;
    pnew->logmarglik = logmarglik;
    pnew->A = new int[lenA];
    for (i = 0; i < lenA; i++)
    {
        pnew->A[i] = A[i];
    }
    p->Next = pnew;
    pnew->Next = pnext;

    // 4. insert the node
    p->Next = pnew;
    pnew->Next = pnext;

    // 5. update counter
    if (counter == nMaxRegs)
    {
        RegressionPointer trash = regression->Next;
        regression->Next = trash->Next;

        delete[] trash->A;
        delete trash;
    }
    else
    {
        counter++;
    }

    return;
}

void DeleteAllRegressions(
    RegressionPointer regressions,
    int &counter)
{

    if (regressions == NULL)
    {
        return;
    }

    RegressionPointer p = regressions->Next;
    RegressionPointer pnext = NULL;

    while (p != NULL)
    {
        // Save next node before deleting current node
        pnext = p->Next;

        // Delete dynamically allocated predictor array
        delete[] p->A;

        p->A = NULL;
        p->Next = NULL;

        // Delete current node
        delete p;

        // Move to next node
        p = pnext;
    }

    // reset head and counter
    regressions->Next = NULL;
    counter = 0;

    return;
}

void SaveRegressions(
    const char *filename,
    RegressionPointer regressions)
{
    SaveRegressions(filename, regressions, NULL);
}

void SaveRegressions(
    const char *filename,
    RegressionPointer regressions,
    const char **columnNames)
{
    if (regressions == NULL)
    {
        printf("Regression list head is NULL.\n");
        return;
    }

    FILE *out = fopen(filename, "w");

    if (out == NULL)
    {
        printf("Cannot open output file [%s]\n", filename);
        exit(1);
    }

    RegressionPointer p = regressions->Next;

    fprintf(out, "logmarglik\tlenA\tpredictors\n");

    while (p != NULL)
    {
        // Print log marginal likelihood and number of predictors
        fprintf(out, "%.10lf\t%d", p->logmarglik, p->lenA);

        // Print predictor indices
        fprintf(out, "\t");
        for (int i = 0; i < p->lenA; i++)
        {
            if (i > 0)
            {
                fprintf(out, ",");
            }
            fprintf(out, "%d", p->A[i]);
        }

        fprintf(out, "\n");

        p = p->Next;
    }

    fclose(out);

    return;
}
