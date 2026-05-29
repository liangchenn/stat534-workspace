#include <stdio.h>

#include "ai_regmodels.h"
#include "matrices.h"

int ai_sameRegression(int lenA, int *A, int lenB, int *B)
{
    if (lenA != lenB)
    {
        return 0;
    }

    for (int i = 0; i < lenA; i++)
    {
        if (A[i] != B[i])
        {
            return 0;
        }
    }

    return 1;
}

int ai_countRegressions(LPAIRegression regressions)
{
    int count = 0;
    LPAIRegression p = regressions->Next;

    while (p != NULL)
    {
        count++;
        p = p->Next;
    }

    return count;
}

void ai_deleteLastRegression(LPAIRegression regressions)
{
    if (regressions->Next == NULL)
    {
        return;
    }

    LPAIRegression prev = regressions;
    LPAIRegression curr = regressions->Next;

    while (curr->Next != NULL)
    {
        prev = curr;
        curr = curr->Next;
    }

    prev->Next = NULL;
    delete[] curr->A;
    delete curr;
}

void ai_addRegression(int nMaxRegs, LPAIRegression regressions, int lenA, int *A, double logmarglikA)
{
    LPAIRegression prev = regressions;
    LPAIRegression curr = regressions->Next;

    while (curr != NULL)
    {
        if (ai_sameRegression(lenA, A, curr->lenA, curr->A))
        {
            return;
        }

        if (curr->logmarglik < logmarglikA)
        {
            break;
        }

        prev = curr;
        curr = curr->Next;
    }

    LPAIRegression node = new AIRegression;
    node->lenA = lenA;
    node->logmarglik = logmarglikA;
    node->A = new int[lenA];

    for (int i = 0; i < lenA; i++)
    {
        node->A[i] = A[i];
    }

    node->Next = curr;
    prev->Next = node;

    if (ai_countRegressions(regressions) > nMaxRegs)
    {
        ai_deleteLastRegression(regressions);
    }
}

void ai_findBestTwoPredictorModels(int nMaxRegs, gsl_matrix *data, int p, LPAIRegression regressions)
{
    int A[2];

    for (int i = 1; i < p; i++)
    {
        A[0] = i + 1;
        ai_addRegression(nMaxRegs, regressions, 1, A, marglik(data, 1, A));
    }

    for (int i = 1; i < p; i++)
    {
        for (int j = i + 1; j < p; j++)
        {
            A[0] = i + 1;
            A[1] = j + 1;
            ai_addRegression(nMaxRegs, regressions, 2, A, marglik(data, 2, A));
        }
    }
}
