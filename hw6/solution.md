# Homework 6 Solution Guide

This file summarizes how to solve each problem in `hw6.pdf` and explains the steps you should follow when writing the code.

## Files in This Folder

- `hw6.pdf`: assignment statement.
- `mybandedmatrix.txt`: the matrix used in Problem 1.
- `example/main.cpp`: example driver for regression model search.
- `example/matrices.cpp`: matrix utilities and marginal likelihood function.
- `example/regmodels.cpp`: linked-list utilities for regression models.

The example code refers to headers such as `matrices.h` and `regmodels.h`, and to the data file `erdata.txt`. Those files are not present in this folder, so the regression output for Problems 2 and 3 cannot be reproduced from the current files alone. The solution below explains the required code changes and what output should be produced once the missing files are available.

## Problem 1: Recursive Determinant

### Goal

Write a recursive determinant function and use it in `main.cpp` to read `mybandedmatrix.txt`, compute the determinant, and print it.

The supplied matrix is `10 x 10`. Its determinant is approximately:

```text
0.5385580062866211
```

### Step 1: Understand the Recursive Formula

For an `n x n` matrix `A`:

1. If `n == 1`, then

   ```text
   det(A) = A[0][0]
   ```

2. If `n == 2`, then

   ```text
   det(A) = A[0][0] * A[1][1] - A[0][1] * A[1][0]
   ```

3. If `n >= 3`, expand along one row. The easiest row is row `0`:

   ```text
   det(A) = sum over j of A[0][j] * (-1)^j * det(minor matrix removing row 0 and column j)
   ```

The assignment uses 1-based notation, so its sign is `(-1)^(i+j)`. In C++ with row `0`, the sign becomes `+ - + - ...`, which is `(-1)^j`.

### Step 2: Build the Minor Matrix

For each column `col` in the first row:

1. Allocate a new `(n - 1) x (n - 1)` matrix.
2. Copy every row except row `0`.
3. Copy every column except `col`.
4. Recursively compute the determinant of that smaller matrix.
5. Free the minor matrix.

### Step 3: Determinant Function

```cpp
double determinant(double** A, int n)
{
  if(n == 1)
  {
    return A[0][0];
  }

  if(n == 2)
  {
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];
  }

  double det = 0.0;

  for(int col = 0; col < n; col++)
  {
    double** minor = allocmatrix(n - 1, n - 1);

    for(int i = 1; i < n; i++)
    {
      int minorCol = 0;
      for(int j = 0; j < n; j++)
      {
        if(j == col)
        {
          continue;
        }

        minor[i - 1][minorCol] = A[i][j];
        minorCol++;
      }
    }

    double sign = (col % 2 == 0) ? 1.0 : -1.0;
    det += sign * A[0][col] * determinant(minor, n - 1);

    freematrix(n - 1, minor);
  }

  return det;
}
```

### Step 4: Main Program for Problem 1

```cpp
int main()
{
  int n = 10;
  char filename[] = "mybandedmatrix.txt";

  double** A = allocmatrix(n, n);
  readmatrix(filename, n, n, A);

  double detA = determinant(A, n);
  printf("determinant = %.16lf\n", detA);

  freematrix(n, A);

  return 0;
}
```

Expected output:

```text
determinant = 0.5385580062866211
```

## Problem 2: Keep Only the Best `nMaxRegs` Regressions

### Goal

Modify `AddRegression` so its signature becomes:

```cpp
void AddRegression(int nMaxRegs,
                   LPRegression regressions,
                   int lenA,
                   int* A,
                   double logmarglikA)
```

The linked list should stay sorted in decreasing order of log marginal likelihood. After inserting a new model, if the list has more than `nMaxRegs` elements, delete the last element. Since the list is sorted from largest to smallest, the last element is the worst model.

### Step 1: Keep the Existing Sorted Insert Logic

The original function already does most of the work:

1. Start at the dummy head node.
2. Walk through the list.
3. If the same regression already exists, return.
4. Stop when the next model has a log marginal likelihood smaller than or equal to the new one.
5. Insert the new model at that position.

This keeps the list sorted from best to worst.

### Step 2: Count the List Length

After inserting the new model, count how many real nodes are in the list:

```cpp
int CountRegressions(LPRegression regressions)
{
  int count = 0;
  LPRegression p = regressions->Next;

  while(p != NULL)
  {
    count++;
    p = p->Next;
  }

  return count;
}
```

### Step 3: Delete the Last Regression if Needed

The example code already provides:

```cpp
void DeleteLastRegression(LPRegression regressions)
```

So after insertion:

```cpp
if(CountRegressions(regressions) > nMaxRegs)
{
  DeleteLastRegression(regressions);
}
```

### Step 4: Modified `AddRegression`

```cpp
void AddRegression(int nMaxRegs,
                   LPRegression regressions,
                   int lenA,
                   int* A,
                   double logmarglikA)
{
  int i;
  LPRegression p = regressions;
  LPRegression pnext = p->Next;

  while(NULL != pnext)
  {
    if(sameregression(lenA, A, pnext->lenA, pnext->A))
    {
      return;
    }

    if(pnext->logmarglikA > logmarglikA)
    {
      p = pnext;
      pnext = p->Next;
    }
    else
    {
      break;
    }
  }

  LPRegression newp = new Regression;
  newp->lenA = lenA;
  newp->logmarglikA = logmarglikA;
  newp->A = new int[lenA];

  for(i = 0; i < lenA; i++)
  {
    newp->A[i] = A[i];
  }

  p->Next = newp;
  newp->Next = pnext;

  if(CountRegressions(regressions) > nMaxRegs)
  {
    DeleteLastRegression(regressions);
  }

  return;
}
```

### Important Header Change

Update the declaration in `regmodels.h` from the old version to:

```cpp
void AddRegression(int nMaxRegs,
                   LPRegression regressions,
                   int lenA,
                   int* A,
                   double logmarglikA);
```

Also declare `CountRegressions` if you define it in `regmodels.cpp`:

```cpp
int CountRegressions(LPRegression regressions);
```

## Problem 3: Best 10 Regressions with At Most Two Predictors

### Goal

Use the modified `AddRegression` to find the 10 regressions with the largest marginal likelihood among all models with one or two predictors.

In the example code:

```cpp
int n = 158;
int p = 51;
```

Column `1` is the response, and predictor columns are `2` through `51`. Because the code uses 1-based column labels for `A`, the valid predictor labels are:

```text
2, 3, ..., 51
```

### Step 1: Set the Maximum Number of Saved Models

```cpp
int nMaxRegs = 10;
```

### Step 2: Add All One-Predictor Models

```cpp
lenA = 1;
for(i = 1; i < p; i++)
{
  A[0] = i + 1;

  AddRegression(nMaxRegs,
                regressions,
                lenA,
                A,
                marglik(n, p, data, lenA, A));
}
```

### Step 3: Add All Two-Predictor Models

Use nested loops with `i < j` so every pair appears once and the predictor vector stays sorted:

```cpp
lenA = 2;
for(i = 1; i < p; i++)
{
  for(j = i + 1; j < p; j++)
  {
    A[0] = i + 1;
    A[1] = j + 1;

    AddRegression(nMaxRegs,
                  regressions,
                  lenA,
                  A,
                  marglik(n, p, data, lenA, A));
  }
}
```

### Step 4: Save the Best 10 Models

```cpp
char outputfilename[] = "best10-regressions1-2.txt";
SaveRegressions(outputfilename, regressions);
```

The output file will contain at most 10 lines. Each line has:

```text
log marginal likelihood    number of predictors    predictor indices
```

For example, the format should look like this:

```text
-123.45678    2    4    17
-124.00342    1    9
```

The actual values require `erdata.txt`, which is not included in this folder.

## Full Main Logic for Problem 3

The important part of `main.cpp` should look like this:

```cpp
int main()
{
  int i, j;

  int n = 158;
  int p = 51;
  int nMaxRegs = 10;

  char datafilename[] = "erdata.txt";
  char outputfilename[] = "best10-regressions1-2.txt";

  double** data = allocmatrix(n, p);
  readmatrix(datafilename, n, p, data);

  LPRegression regressions = new Regression;
  regressions->Next = NULL;

  int A[p - 1];
  int lenA;

  lenA = 1;
  for(i = 1; i < p; i++)
  {
    A[0] = i + 1;

    AddRegression(nMaxRegs,
                  regressions,
                  lenA,
                  A,
                  marglik(n, p, data, lenA, A));
  }

  lenA = 2;
  for(i = 1; i < p; i++)
  {
    for(j = i + 1; j < p; j++)
    {
      A[0] = i + 1;
      A[1] = j + 1;

      AddRegression(nMaxRegs,
                    regressions,
                    lenA,
                    A,
                    marglik(n, p, data, lenA, A));
    }
  }

  SaveRegressions(outputfilename, regressions);

  DeleteAllRegressions(regressions);
  freematrix(n, data);
  delete regressions;
  regressions = NULL;

  return 0;
}
```

## How to Think Through the Assignment

1. Problem 1 is about recursion. Each recursive call solves the same determinant problem on a smaller matrix. The base cases stop the recursion at `1 x 1` or `2 x 2`.

2. Problem 2 is about maintaining an ordered linked list. You do not need to sort the whole list after every insertion. Insert the new regression directly into the correct position, then remove the tail if the list is too long.

3. Problem 3 is about enumeration. Generate every model with one predictor, then every model with two predictors. Call `AddRegression(10, ...)` each time. Because `AddRegression` keeps only the best 10, you never need to store all candidate models at once.

4. The key invariant is: after every call to `AddRegression`, the linked list is sorted from largest to smallest marginal likelihood and contains no more than `nMaxRegs` models.

