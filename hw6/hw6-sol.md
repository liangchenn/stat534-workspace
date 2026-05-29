# STAT534 Hw6

## Problem 1: Recursive Determinant

### Part A: My Solution

- I use the Laplace decomposition provided by the problem set, and expand at the first row (i.e. `i == 0`).
- I use recursion to calculate the determinant with only one base 1-dim case, since 2-dim case is covered by the algorithm.
- Here I accelerate the recursion by skipping `aij = 0` cases, since it won't affect the determinant.
- I use GSL lib to manage the matrix mem and operation, use  `getSubmatrix` to create each sub matrix and frees it after use.
- Output: `The determinant of mybandedmatrix.txt is 0.5385580063.`


```cpp
double getDeterminant(gsl_matrix *matrix)
{
    double det = 0.0;
    int n_rows = matrix->size1;
    int n_cols = matrix->size2;

    if (n_rows != n_cols)
    {
        printf("Error: matrix must be square, currently got (%d, %d).\n", n_rows, n_cols);
        exit(1);
    }

    if (n_rows == 1)
    {
        return gsl_matrix_get(matrix, 0, 0);
    }

    const int i = 0;

    for (int j = 0; j < n_cols; j++)
    {
        double aij = gsl_matrix_get(matrix, i, j);
        if (aij == 0)
        {
            continue;
        }

        gsl_matrix *submatrix = getSubmatrix(matrix, i, j);
        double subdet = getDeterminant(submatrix);
        det += ((j % 2 == 0 ? 1.0 : -1.0) * aij * subdet);
        gsl_matrix_free(submatrix);
    }

    return det;
}
```


### Part B: AI Solution

#### Prompt

Read Problem 1 from `hw6.pdf`. Write C++ code recursively calculate the determinant of an `n x n` matrix. You could refer to the example matrices.cpp codes to performa matrix operation. Use base cases for `n = 1` and `n = 2`, and for `n >= 3` use Laplace expansion with submatrices formed by removing one row and one column. Put the AI code in `ai_matrices.h` and `ai_matrices.cpp`.

#### Code Generated

```cpp
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
```

### Part C: Comparison

- **Overall**
  - Both correct solutions use recursive determinant calculation.
  - AI solution is a little more complete because it has both `n = 1` and `n = 2` base cases, but in my opinion redundant.

- **Correctness**
  - Both expand along the first row and use the same sign pattern.
  - Both give the correct determinant = 0.5385580063.

- **Data Structure**
  - Both use `gsl_matrix*` for handling matrix data and operations.

- **Readability**
  - My code has clearer comments for the homework steps.
  - My code is shorter and has a direct `n = 2` formula.

- **Algorithm**
  - Both are Laplace expansion algorithms.
  - Both skip zero entries, which is faster for the banded matrix.



## Problem 2: AddRegression

### Part A: My Solution

- I modify the variable and struct name to be more informative.
- I modify the `AddRegression` function to accept a global variable, `counter`, to keep track on the number of the models in the linked list.
- I modify the linked list order to be ascending instead of descending, since I don't want to traverse all models to delete a node. Also, when the list is full, it's faster to reject a smaller model to be inserted.
- I also remove the same-model-check logic from the function, since I think it could be easily handled outside the function with well-defined for loop in this case (See problem 3).

```cpp
// my solution
void AddRegression(
    int nMaxRegs,
    RegressionPointer regression,
    int lenA,
    int *A,
    double logmarglik,
    int &counter)
{
    if (counter >= nMaxRegs &&
        regression->Next != NULL &&
        logmarglik <= regression->Next->logmarglik)
    {
        // skip if full and loglik is smaller than the smallest in the list
        return;
    }

    // create 2 pointers to keep track curr & next node
    RegressionPointer p = regression;
    RegressionPointer pnext = p->Next;
    // find insert position
    while (pnext != NULL)
    {
        if (logmarglik > pnext->logmarglik)
        {
            // move forward 2 pointers
            p = pnext;
            pnext = p->Next;
        }
        else
        {
            break;
        }
    }

    // create new model node, insert to the position
    RegressionPointer pnew = new Regression;
    pnew->lenA = lenA;
    pnew->logmarglik = logmarglik;
    pnew->A = new int[lenA];
    for (int i = 0; i < lenA; i++)
    {
        pnew->A[i] = A[i];
    }
    p->Next = pnew;
    pnew->Next = pnext;

    // delete & update logic
    if (counter == nMaxRegs)
    {
        RegressionPointer trash = regression->Next;
        regression->Next = trash->Next;
        delete[] trash->A;
        delete trash;
    }
    else
    {   
        // update counter if new node inserted
        counter++;
    }
}
```

Original Code:
- descending order linked list
    - bad for delete smallest loglik model imo.
- variable and struct alias a bit ambiguous.

```cpp
void AddRegression(LPRegression regressions,int lenA,int* A,double logmarglikA)
{
  int i;
  LPRegression p = regressions;
  LPRegression pnext = p->Next;

  while(NULL!=pnext)
  {
     //return if we have previously found this regression
     if(sameregression(lenA,A,pnext->lenA,pnext->A))
     {
        return;
     }

     //go to the next element in the list if the current
     //regression has a larger log marginal likelihood than
     //the new regression A
     if(pnext->logmarglikA>logmarglikA)
     {
        p = pnext;
        pnext = p->Next;
     }
     else //otherwise stop; this is where we insert the new regression
     {
        break;
     }
  }

  //create a new element of the list
  LPRegression newp = new Regression;
  newp->lenA = lenA;
  newp->logmarglikA = logmarglikA;
  newp->A = new int[lenA];
  
  //copy the predictors
  for(i=0;i<lenA;i++)
  {
    newp->A[i] = A[i];
  }

  //insert the new element in the list
  p->Next = newp;
  newp->Next = pnext;

  printf("inserted [%d]\n",A[0]);

  return;
}
```


### Part B: AI Solution

#### Prompt

Read Problem 2 from `hw6.pdf`. Transform `AddRegression` so it receives `nMaxRegs`, a linked list of regressions, the predictor set `A`, and `logmarglikA`. The function should keep only the `nMaxRegs` regressions with largest marginal likelihood. Put the AI code in `ai_regmodels.h` and `ai_regmodels.cpp`.

#### Code Generated

```cpp
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
```

### Part C: Comparison

- **Overall**
  - Both solutions keep only a fixed number of best regressions.
  - My solution is more efficient regarding time complexity overall.
  - For top-k models, my solution uses a counter, while AI counts and deletes the last node.

- **Correctness**
  - Both  solutions give the correct top 10 models for this program.
  - AI solution follows the assignment wording more directly by keeping the list sorted best-to-worst.
  - AI also checks duplicate regressions explicitly, but mine left this to outer logic. (But sameRegression function is really slow; and it also cannot check the past model that is deleted, which is a bit inefficient.)

- **Data Structure**
  - Both use a linked list and dynamically allocated predictor arrays.
  - My solution uses `RegressionPointer`; AI uses a separate `LPAIRegression` to avoid conflict.

- **Readability**
  - My solution is compact, and the main logic is self-contained, I think it's more readable.
  - AI solution is more modular because counting, duplicate check, and deletion are separate helper functions.

- **Algorithm**
  - Both insert one model at a time in sorted position.
  - My solution avoids recounting the list, AI's need to check everytime, and could only check 10 models at most at one time. It's possible to re-count the previous (but deleted) models.
  - AI solution is simpler but may be slower because it counts the list after insertion.
  - Complexity:
    - Deletion:
        - Mine is O(1); AI's is O(k)
    - Top-K logic:
        - Mine is handled by global var., hence O(1); AI's check every time traversing the list, O(k)
    - Insertion:
        - Mine: O(k); AI's: O(k)


## Problem 3: Best 10 Regressions

### Part A: My Solution

- Enumerates all one-predictor models.
- Enumerates all two-predictor models with `i < j`, so there are no repeated pairs.
- Uses `AddRegression` to keep only the best 10 models instead of storing all models.


```cpp
// 0. constants
    const int nMaxRegressions = 10;

    // 1. load data
    int i, j;

    int n = 158;                        // sample size
    int p = 51;                         // number of variables
    char datafilename[] = "erdata.txt"; // name of the data file
    char outputfilename[] = "top-10-regression-models.txt";
    char columnNameStorage[51][16];
    const char *columnNames[51];

    snprintf(columnNameStorage[0], sizeof(columnNameStorage[0]), "Y");
    columnNames[0] = columnNameStorage[0];
    for (i = 1; i < p; i++)
    {
        snprintf(columnNameStorage[i], sizeof(columnNameStorage[i]), "X%d", i);
        columnNames[i] = columnNameStorage[i];
    }

    // allocate the data matrix
    gsl_matrix *data = gsl_matrix_alloc(n, p);

    // read the data
    FILE *datafile = fopen(datafilename, "r");
    // readmatrix(datafilename, n, p, data);
    gsl_matrix_fscanf(datafile, data);
    fclose(datafile);

    // 2. Add regression model for sorted, non-dup. regressors

    // (0) max counter
    int counter = 0;

    // (1) create linked list
    RegressionPointer regressions = new Regression;

    regressions->lenA = 0;
    regressions->logmarglik = 0.0;
    regressions->A = NULL;
    regressions->Next = NULL;

    // (2) single predictor case
    int lenA = 1;
    int A[2];
    for (i = 1; i < p; i++) // first column is response
    {
        A[0] = i + 1;
        double logmarglik = marglik(data, lenA, A);
        AddRegression(nMaxRegressions, regressions, lenA, A, logmarglik, counter);
    }

    // (3) two-predictor case
    lenA = 2;

    for (i = 1; i < p; i++)
    {
        for (j = i + 1; j < p; j++)
        {
            // create non-duplicated index pair by only including {i, j > i}
            A[0] = i + 1;
            A[1] = j + 1;

            double logmarglik = marglik(data, lenA, A);
            AddRegression(nMaxRegressions, regressions, lenA, A, logmarglik, counter);
        }
    }

    // (4) output result
    SaveRegressions(outputfilename, regressions, columnNames);

    // (5) release mem
    DeleteAllRegressions(regressions, counter);
    gsl_matrix_free(data);
```

Output:

```text
logmarglik       lenA    predictors
-69.5779357816   2      4,24
-68.3023725325   2      2,33
-67.9188394199   2      2,8
-67.3221733280   2      2,12
-66.2808759392   2      2,34
-65.6531121432   2      2,24
-65.3939518081   2      2,15
-64.1045542949   2      2,43
-61.9074377340   2      2,10
-59.3162035297   2      2,11
```

### Part B: AI Solution

#### Prompt

Read Problem 3 from `hw6.pdf`. Use the modified `AddRegression` from Problem 2 to find the best 10 regressions with at most two predictors. The data has 158 rows and 51 columns; column 1 is the response and columns 2 to 51 are predictors. Generate all one-predictor and two-predictor models, compute marginal likelihood by `marglik`, and retain the 10 largest values. Put the AI code in `ai_regmodels.h` and `ai_regmodels.cpp`.

#### Code Generated

```cpp
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
```

### Part C: Comparison

- **Overall**
  - Both solutions enumerate exactly the required model space.
  - My solution is faster than AI's due to the time complexity difference.

- **Correctness**
  - Both use one-based predictor labels through `i + 1`.
  - Both use `i < j` for two-predictor models.
  - Both rely on `AddRegression` to retain the largest marginal likelihoods.

- **Data Structure**
  - Both use a fixed local array `A[2]` for candidate predictors.
  - Both use a linked list to store the current top 10 models.

- **Readability**
  - Both is concise regarding the main logic.

- **Algorithm**
  - The main algo. in problem 3 is how we construct the non-duplicated index pairs. Both solutions choose to use for-loop with i < j condition.
  - The AI solution is slower due to slow algo in problem 2's function.
