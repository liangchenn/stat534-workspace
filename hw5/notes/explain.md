# matrices.cpp / matrices.h 中文說明

這兩個檔案合起來是一個「矩陣工具庫」：

- `matrices.h` 是 header file，負責宣告可被其他 `.cpp` 檔使用的函式。
- `matrices.cpp` 是 implementation file，負責實作這些函式。
- `main.cpp` 只要 `#include "matrices.h"`，就可以呼叫 `allocmatrix`、`matrixproduct`、`inverse`、`logdet` 等函式。

這份程式使用的是偏 C 風格的 C++ 寫法：矩陣用 `double**` 表示，檔案 I/O 用 `FILE*`、`fopen`、`fscanf`、`fprintf`，錯誤處理用 `printf` / `fprintf` 加上 `exit(1)`。同時，它也用了 C++ 的 `new[]` / `delete[]` 來配置和釋放二維矩陣。

## 整體設計

矩陣在這份程式中的主要表示方式是：

```cpp
double** m;
```

可以把它理解成「指向很多 row 的指標」：

```text
m
|
|-- m[0] -> 第 0 row: double[p]
|-- m[1] -> 第 1 row: double[p]
|-- m[2] -> 第 2 row: double[p]
...
```

所以如果矩陣有 `n` 列、`p` 欄：

```cpp
m[i][j]
```

代表第 `i` 列、第 `j` 欄的元素。C/C++ 陣列從 0 開始，所以 `i` 的合法範圍是 `0` 到 `n - 1`，`j` 的合法範圍是 `0` 到 `p - 1`。

這種設計的特色：

- 直覺，因為可以用 `m[i][j]` 存取矩陣元素。
- 每一列分別配置記憶體，不一定連續。
- 必須自己手動釋放記憶體，否則會 memory leak。
- 和 LAPACK/LAPACKE 這類線性代數 library 溝通時，常常需要轉成一維陣列。

## matrices.h 的作用

`matrices.h` 是標頭檔，內容大致分成兩類：

```cpp
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <lapacke.h>
```

這些是函式會用到的標準 library 或外部 library。

- `<stdio.h>`：提供 `FILE*`、`fopen`、`fscanf`、`fprintf`、`printf`。
- `<math.h>`：提供 `log`、`fabs` 等數學函式。
- `<string.h>`：提供 `memset`。
- `<stdlib.h>`：提供 `malloc`、`free`、`exit`。
- `<lapacke.h>`：提供 LAPACKE 線性代數函式。

接著是函式宣告：

```cpp
double ** allocmatrix(int n,int p);
void freematrix(int n,double** m);
void copymatrix(int n,int p,double** source,double** dest);
void readmatrix(char* filename,int n,int p,double* m[]);
void printmatrix(char* filename,int n,int p,double** m);
double** transposematrix(int n,int p,double** m);
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m);
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m);
void inverse(int p,double** m);
double logdet(int p,double** m);
```

這些宣告的意思是：「別的檔案可以知道這些函式存在，以及要怎麼呼叫它們。」例如 `main.cpp` 裡面只要：

```cpp
#include "matrices.h"
```

編譯器就知道 `allocmatrix` 是一個接受兩個 `int`、回傳 `double**` 的函式。

### Header file 的 convention

一般 C/C++ 專案常見的 header file 會有 include guard：

```cpp
#ifndef MATRICES_H
#define MATRICES_H

// declarations

#endif
```

或使用：

```cpp
#pragma once
```

這份 `matrices.h` 目前沒有 include guard。小作業通常還可以，但大型專案中可能造成重複 include 的問題。

另外這份 header 手動宣告了：

```cpp
int LAPACKE_dgeev(...);
```

但因為已經 `#include <lapacke.h>`，通常不需要自己重新宣告 LAPACKE 函式。重複宣告如果 signature 不完全一致，可能造成編譯或連結問題。

## matrices.cpp 的作用

`matrices.cpp` 是實作檔，第一行：

```cpp
#include "matrices.h"
```

表示它引入自己的 header，確保 implementation 和 declaration 能對上。

第二行：

```cpp
#include <lapacke.h>
```

再次引入 LAPACKE。因為 `matrices.h` 已經 include 了 `<lapacke.h>`，這行不一定必要，但通常不會造成問題。

## 函式逐一說明

### `allocmatrix`

```cpp
double ** allocmatrix(int n,int p)
{
    int i;
    double** m;

    m = new double*[n];
    for(i=0;i<n;i++)
    {
        m[i] = new double[p];
        memset(m[i],0,p*sizeof(double));
    }
    return(m);
}
```

功能：配置一個 `n x p` 的矩陣，並把所有元素初始化成 0。

重要語法：

```cpp
double** m;
```

`m` 是 pointer to pointer，也就是「指向 `double*` 的指標」。這裡用來表示二維矩陣。

```cpp
m = new double*[n];
```

配置一個長度為 `n` 的陣列，每個元素都是 `double*`。也就是先配置 row pointers。

```cpp
m[i] = new double[p];
```

對每一列配置 `p` 個 `double`。

```cpp
memset(m[i],0,p*sizeof(double));
```

把第 `i` 列的記憶體全部設成 0。`sizeof(double)` 是一個 double 佔用的 byte 數，`p * sizeof(double)` 就是整列的 byte 數。

呼叫方式：

```cpp
double** A = allocmatrix(3, 4);
```

這會得到一個 3 列 4 欄的矩陣，所有值一開始都是 0。

### `freematrix`

```cpp
void freematrix(int n,double** m)
{
    int i;

    for(i=0;i<n;i++)
    {
        delete[] m[i]; m[i] = NULL;
    }
    delete[] m; m = NULL;
    return;
}
```

功能：釋放由 `allocmatrix` 配置的矩陣。

因為 `allocmatrix` 是先配置 row pointer array，再分別配置每一列，所以釋放時也要反過來：

1. 先對每一列 `delete[] m[i]`
2. 再對 row pointer array `delete[] m`

重要 convention：

- `new[]` 要搭配 `delete[]`。
- `malloc` 要搭配 `free`。
- 不要用 `free` 釋放 `new[]` 出來的記憶體。
- 不要用 `delete[]` 釋放 `malloc` 出來的記憶體。

```cpp
m[i] = NULL;
m = NULL;
```

這是在函式內把指標設成空指標，避免 dangling pointer。不過要注意，`m = NULL` 只改到函式內部的 local copy，呼叫端的指標不會變成 `NULL`。

### `copymatrix`

```cpp
void copymatrix(int n,int p,double** source,double** dest)
{
    int i,j;

    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            dest[i][j] = source[i][j];
        }
    }
    return;
}
```

功能：把 `source` 矩陣複製到 `dest` 矩陣。

設計上它應該要複製 `n x p` 個元素，所以外層跑列 `i < n`，內層應該跑欄 `j < p`。

目前程式寫成：

```cpp
for(j=0;j<n;j++)
```

這是一個 bug。若矩陣是方陣，例如 `n == p`，暫時看不出問題；但如果是非方陣，例如 `2 x 5`，它只會複製前 2 欄。如果是 `5 x 2`，它還可能越界存取。

比較合理的版本：

```cpp
for(j=0;j<p;j++)
```

### `readmatrix`

```cpp
void readmatrix(char* filename,int n,int p,double* m[])
{
    int i,j;
    double s;
    FILE* in = fopen(filename,"r");

    if(NULL==in)
    {
        printf("Cannot open input file [%s]\n",filename);
        exit(1);
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<p;j++)
        {
            fscanf(in,"%lf",&s);
            m[i][j] = s;
        }
    }
    fclose(in);
    return;
}
```

功能：從文字檔讀入 `n x p` 個 double，存入矩陣 `m`。

重要語法：

```cpp
FILE* in = fopen(filename,"r");
```

開啟檔案。`"r"` 表示 read mode。

```cpp
if(NULL==in)
```

檢查檔案是否成功開啟。`NULL == in` 和 `in == NULL` 意思一樣。有些 C 程式設計師喜歡把常數放左邊，避免不小心寫成 assignment。

```cpp
fscanf(in,"%lf",&s);
```

從檔案讀一個 double 到 `s`。`%lf` 是 `double` 的格式，`&s` 表示傳入 `s` 的地址，讓 `fscanf` 可以改它。

```cpp
m[i][j] = s;
```

把讀到的值放進矩陣。

注意：這裡沒有檢查 `fscanf` 是否真的成功讀到數字。更穩健的寫法會檢查回傳值是否為 1。

### `printmatrix`

```cpp
void printmatrix(char* filename,int n,int p,double** m)
{
    int i,j;
    double s;
    FILE* out = fopen(filename,"w");

    if(NULL==out)
    {
        printf("Cannot open output file [%s]\n",filename);
        exit(1);
    }
    for(i=0;i<n;i++)
    {
        fprintf(out,"%.3lf",m[i][0]);
        for(j=1;j<p;j++)
        {
            fprintf(out,"\t%.3lf",m[i][j]);
        }
        fprintf(out,"\n");
    }
    fclose(out);
    return;
}
```

功能：把矩陣輸出到文字檔。

重要語法：

```cpp
FILE* out = fopen(filename,"w");
```

用 write mode 開啟檔案。如果檔案已存在，內容會被覆蓋。

```cpp
fprintf(out,"%.3lf",m[i][0]);
```

把 double 印到檔案，保留小數點後 3 位。

```cpp
fprintf(out,"\t%.3lf",m[i][j]);
```

`\t` 是 tab，用來分隔欄位。

這個函式中：

```cpp
double s;
```

宣告了但沒有使用，可以刪掉。

### `transposematrix`

```cpp
double** transposematrix(int n,int p,double** m)
{
    int i,j;

    double** tm = allocmatrix(p,n);

    for(i=0;i<p;i++)
    {
        for(j=0;j<n;j++)
        {
            tm[i][j] = m[j][i];
        }
    }

    return(tm);
}
```

功能：建立並回傳矩陣轉置。

如果原矩陣 `m` 是 `n x p`，轉置後 `tm` 是 `p x n`。

數學上：

```text
tm[i][j] = m[j][i]
```

例如：

```text
m = [1 2 3
     4 5 6]

transpose(m) = [1 4
                2 5
                3 6]
```

重要設計：

```cpp
double** tm = allocmatrix(p,n);
return(tm);
```

這個函式內部配置了新矩陣並回傳，所以呼叫者之後必須呼叫：

```cpp
freematrix(p, tm);
```

否則會 memory leak。

### `dotmatrixproduct`

```cpp
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m)
{
    int i,j;

    for(i=0;i<n;i++)
    {
        for(j=0;j<p;j++)
        {
            m[i][j] = m1[i][j]*m2[i][j];
        }
    }

    return;
}
```

功能：做 element-wise product，也就是每個位置各自相乘。

數學上：

```text
m[i][j] = m1[i][j] * m2[i][j]
```

這不是一般線性代數裡的 matrix multiplication，而是 Hadamard product。

例如：

```text
[1 2]   [10 20]   [10 40]
[3 4] * [30 40] = [90 160]
```

設計 pattern：

```cpp
void dotmatrixproduct(..., double** m)
```

它不回傳新矩陣，而是把結果寫進呼叫者傳入的 `m`。所以使用前要先配置好：

```cpp
double** C = allocmatrix(n, p);
dotmatrixproduct(n, p, A, B, C);
```

### `matrixproduct`

```cpp
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m)
{
    int i,j,k;
    double s;

    for(i=0;i<n;i++)
    {
        for(k=0;k<l;k++)
        {
            s = 0;
            for(j=0;j<p;j++)
            {
                s += m1[i][j]*m2[j][k];
            }
            m[i][k] = s;
        }
    }
    return;
}
```

功能：做一般矩陣乘法。

輸入：

- `m1` 是 `n x p`
- `m2` 是 `p x l`
- `m` 是結果矩陣，大小是 `n x l`

數學上：

```text
m[i][k] = sum over j of m1[i][j] * m2[j][k]
```

三層迴圈的角色：

- `i`：結果矩陣的 row。
- `k`：結果矩陣的 column。
- `j`：內積中間的 summation index。

例子：

```cpp
double** C = allocmatrix(n, l);
matrixproduct(n, p, l, A, B, C);
```

這個函式假設使用者傳入的矩陣尺寸正確。它不會自己檢查維度。

### `set_mat_identity`

```cpp
void set_mat_identity(int p, double *A)
{
    int i;

    for(i = 0; i < p * p; i++) A[i] = 0;
    for(i = 0; i < p; i++) A[i * p + i] = 1;
    return;
}
```

功能：把一個一維陣列 `A` 設成 `p x p` identity matrix。

identity matrix 長這樣：

```text
[1 0 0
 0 1 0
 0 0 1]
```

重要語法：

```cpp
double *A
```

這裡不是 `double**`，而是一維陣列。LAPACK 常用一維陣列表示矩陣。

```cpp
A[i * p + i] = 1;
```

這是在一維陣列裡設定對角線元素。

注意：這個函式沒有出現在 `matrices.h`，所以目前比較像 `matrices.cpp` 內部使用的 helper function。

### `inverse`

```cpp
void inverse(int p, double** m) {
    int i, j, k, info;
    double* m_copy = (double*)malloc((p * p) * sizeof(double));
    double* m_inv = (double*)malloc((p * p) * sizeof(double));

    for (j = 0; j < p; j++) {
        for (i = 0; i < p; i++) {
            m_copy[j * p + i] = m[i][j];
        }
    }

    set_mat_identity(p, m_inv);

    info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', p, p, m_copy, p, m_inv, p);

    if (info != 0) {
        fprintf(stderr, "Something was wrong with LAPACKE_dposv [%d]\n", info);
        exit(1);
    }

    for (j = 0; j < p; j++) {
        for (i = 0; i < p; i++) {
            m[i][j] = m_inv[j * p + i];
        }
    }

    free(m_copy);
    free(m_inv);
}
```

功能：計算一個 `p x p` 對稱正定矩陣的反矩陣，並直接覆蓋原本的 `m`。

也就是呼叫：

```cpp
inverse(p, A);
```

之後，`A` 本身會變成 `A^{-1}`。這種設計叫做 in-place modification。

#### 為什麼要轉成一維陣列？

原本矩陣是：

```cpp
double** m;
```

但 LAPACKE 通常吃的是：

```cpp
double* a;
```

也就是一段連續記憶體。所以程式配置：

```cpp
double* m_copy = (double*)malloc((p * p) * sizeof(double));
double* m_inv = (double*)malloc((p * p) * sizeof(double));
```

`m_copy` 存原矩陣，`m_inv` 一開始存 identity matrix，之後會變成反矩陣。

#### Column-major order

這段：

```cpp
m_copy[j * p + i] = m[i][j];
```

是在把 `double**` 矩陣轉成 LAPACK 常用的 column-major order。

C/C++ 一般比較常用 row-major：

```text
index = row * number_of_columns + column
```

但 Fortran / LAPACK 傳統上用 column-major：

```text
index = column * number_of_rows + row
```

因為函式呼叫中使用：

```cpp
LAPACK_COL_MAJOR
```

所以資料也要用 column-major 排列。

#### LAPACKE_dposv

```cpp
info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', p, p, m_copy, p, m_inv, p);
```

`dposv` 是用來解對稱正定矩陣的線性系統：

```text
A X = B
```

如果把 `B` 設成 identity matrix `I`，那麼：

```text
A X = I
```

解出來的 `X` 就是：

```text
A^{-1}
```

所以這個函式不是直接呼叫「inverse routine」，而是透過解線性方程組求反矩陣。這是常見做法。

參數 `'U'` 表示使用矩陣的 upper triangular part。因為對稱矩陣上下三角互相對稱，LAPACK 只需要讀一半。

#### 錯誤處理

```cpp
if (info != 0) {
    fprintf(stderr, "Something was wrong with LAPACKE_dposv [%d]\n", info);
    exit(1);
}
```

`info == 0` 表示成功。非 0 表示失敗，例如矩陣不是正定矩陣或參數錯誤。

#### 注意事項

- `inverse` 只適合對稱正定矩陣。
- 它會覆蓋原矩陣。
- 如果還需要原矩陣，呼叫前應該先 copy。
- `int k` 被宣告但沒有實際使用，可以刪除。
- 使用 `malloc` 配置，所以最後用 `free` 釋放，這點是正確搭配。

### `logdet`

```cpp
double logdet(int p, double** m) {
    int i, j, info;
    double* a = (double*)malloc(p * p * sizeof(double));
    double* wr = (double*)malloc(p * sizeof(double));
    double* wi = (double*)malloc(p * sizeof(double));

    for (i = 0; i < p; i++) {
        for (j = 0; j < p; j++) {
            a[j * p + i] = m[i][j];
        }
    }

    info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', p, a, p, wr, wi, NULL, 1, NULL, 1);

    if (info != 0) {
        printf("Error in eigenvalue computation [info = %d]\n", info);
        exit(1);
    }

    double logdet = 0.0;
    for (i = 0; i < p; i++) {
        logdet += log(fabs(wr[i]));
    }

    free(a);
    free(wr);
    free(wi);

    return logdet;
}
```

功能：計算矩陣 determinant 的 log，也就是：

```text
log(det(A))
```

它的做法是先計算 eigenvalues，然後利用：

```text
det(A) = product of eigenvalues
log(det(A)) = sum of log(eigenvalues)
```

程式中：

```cpp
wr
```

存 eigenvalue 的 real part。

```cpp
wi
```

存 eigenvalue 的 imaginary part。

```cpp
LAPACKE_dgeev(..., 'N', 'N', ...)
```

`dgeev` 是一般矩陣的 eigenvalue routine。兩個 `'N'` 表示不要計算 left eigenvectors，也不要計算 right eigenvectors，只要 eigenvalues。

```cpp
logdet += log(fabs(wr[i]));
```

這行把每個 real eigenvalue 取絕對值後取 log，再加總。

#### 注意事項

如果矩陣是對稱正定矩陣，所有 eigenvalues 都應該是正實數，所以比較合理的數學式是：

```cpp
logdet += log(wr[i]);
```

目前使用 `fabs(wr[i])` 可以避免負值造成 `log` 無法處理，但也可能掩蓋矩陣不是正定或 determinant 為負的問題。

另外，這個函式配置了 `wi`，但最後完全沒有檢查 imaginary part 是否接近 0。若矩陣不是對稱矩陣，eigenvalues 可能是複數，這時只看 `wr` 會不完整。

對稱正定矩陣更常見也更穩定的 `logdet` 做法是 Cholesky decomposition：

```text
A = L L^T
logdet(A) = 2 * sum(log(diagonal entries of L))
```

但這份程式目前是用 eigenvalue 方法。

## 常見 C/C++ 語法整理

### `#include`

```cpp
#include "matrices.h"
#include <stdio.h>
```

- `#include "..."` 通常用來 include 自己專案裡的檔案。
- `#include <...>` 通常用來 include 系統或 library 的 header。

### Function declaration vs definition

在 `matrices.h`：

```cpp
double ** allocmatrix(int n,int p);
```

這是 declaration，只告訴編譯器函式長什麼樣子。

在 `matrices.cpp`：

```cpp
double ** allocmatrix(int n,int p)
{
    ...
}
```

這是 definition，真正寫出函式要做什麼。

### Return type

```cpp
double** allocmatrix(...)
```

回傳一個 `double**`。

```cpp
void freematrix(...)
```

`void` 表示不回傳東西。

```cpp
double logdet(...)
```

回傳一個 `double`。

### Pointer

```cpp
double* a;
```

`a` 是指向 double 的指標，常用來代表一維 double array。

```cpp
double** m;
```

`m` 是指向 `double*` 的指標，這份程式用它代表二維矩陣。

### Array indexing

```cpp
m[i][j]
```

第 `i` 列第 `j` 欄。

```cpp
a[j * p + i]
```

用一維陣列表示二維矩陣的某個元素。這裡是 column-major order。

### `for` loop

```cpp
for(i=0;i<n;i++)
{
    ...
}
```

三個部分：

```text
初始化：i = 0
條件：i < n
每次迴圈後更新：i++
```

`i++` 表示 `i = i + 1`。

### `if`

```cpp
if(NULL==in)
{
    ...
}
```

如果條件成立，就執行大括號內的程式。

### `return`

```cpp
return(m);
```

回傳 `m`。

```cpp
return;
```

在 `void` 函式裡提早結束函式。放在函式最後時其實可以省略。

### Memory allocation

這份程式同時用了 C++ 和 C 的配置方式。

C++：

```cpp
new double[p];
delete[] m[i];
```

C：

```cpp
malloc(p * p * sizeof(double));
free(a);
```

兩套不要混用。

### `char* filename`

```cpp
void readmatrix(char* filename, ...)
```

`char*` 是 C-style string，也就是字元陣列的指標。現代 C++ 常見寫法會用：

```cpp
const char* filename
```

或：

```cpp
std::string filename
```

因為 `readmatrix` 不會修改 filename，所以 `const char*` 會更準確。

## Function design pattern

這份程式裡有幾種函式設計模式。

### 1. Allocate and return

```cpp
double** allocmatrix(int n,int p);
double** transposematrix(int n,int p,double** m);
```

這類函式會配置新記憶體，然後把 pointer 回傳。呼叫者要負責釋放。

使用 pattern：

```cpp
double** A = allocmatrix(n, p);
...
freematrix(n, A);
```

### 2. Write result into output argument

```cpp
void matrixproduct(..., double** m);
void dotmatrixproduct(..., double** m);
void copymatrix(..., double** dest);
```

這類函式不回傳結果，而是把結果寫進最後傳入的矩陣。

使用 pattern：

```cpp
double** C = allocmatrix(n, l);
matrixproduct(n, p, l, A, B, C);
```

好處是呼叫者可以控制記憶體配置。缺點是呼叫者必須確保 `C` 的尺寸正確。

### 3. Modify in place

```cpp
void inverse(int p,double** m);
```

這類函式會直接修改傳入的矩陣。

使用 pattern：

```cpp
double** Ainv = allocmatrix(p, p);
copymatrix(p, p, A, Ainv);
inverse(p, Ainv);
```

如果直接：

```cpp
inverse(p, A);
```

那原本的 `A` 就不見了，變成 `A^{-1}`。

### 4. Compute and return scalar

```cpp
double logdet(int p,double** m);
```

這類函式不配置回傳矩陣，而是回傳一個數值。

使用 pattern：

```cpp
double ld = logdet(p, A);
```

## 這份程式的功能總結

`matrices.cpp` / `matrices.h` 提供以下功能：

- 配置矩陣：`allocmatrix`
- 釋放矩陣：`freematrix`
- 複製矩陣：`copymatrix`
- 從檔案讀矩陣：`readmatrix`
- 把矩陣寫到檔案：`printmatrix`
- 矩陣轉置：`transposematrix`
- element-wise product：`dotmatrixproduct`
- 一般矩陣乘法：`matrixproduct`
- 對稱正定矩陣反矩陣：`inverse`
- determinant 的 log：`logdet`

## main.cpp 如何使用它

`main.cpp` 中目前有一段被註解掉的測試流程：

1. 建立一個 `2 x 2` 矩陣 `A`。
2. 設定：

```text
A = [4 1
     1 3]
```

3. 呼叫 `logdet(p, A)`。
4. 先 copy `A` 到 `Ainv`，因為 `inverse` 會改掉矩陣。
5. 呼叫 `inverse(p, Ainv)`。
6. 用 `matrixproduct(A, Ainv)` 檢查是否接近 identity matrix。
7. 呼叫 `freematrix` 釋放所有矩陣。

這是一個合理的測試：如果 `A * inverse(A)` 接近 identity matrix，代表 inverse 的結果大致正確。

## makefile 和編譯方式

`makefile` 中：

```makefile
matrices.o: matrices.cpp matrices.h
	gcc -g -c matrices.cpp -o matrices.o ${MATRICES_INCLUDE}
```

這行把 `matrices.cpp` 編譯成 object file。

```makefile
main.o: main.cpp matrices.h
	gcc -g -c main.cpp -o main.o ${MATRICES_INCLUDE}
```

這行把 `main.cpp` 編譯成 object file。

```makefile
matrices: main.o matrices.o
	gcc main.o matrices.o -o matrices ${MATRICES_LIB} ${MATRICES_LSTDFLG}
```

這行把兩個 object files link 成可執行檔 `matrices`。

```makefile
MATRICES_LSTDFLG = -lstdc++ -lopenblas -lm
```

代表 link 時需要：

- `-lstdc++`：C++ standard library。
- `-lopenblas`：OpenBLAS / LAPACK library。
- `-lm`：math library。

雖然檔案是 `.cpp`，makefile 使用的是 `gcc`，然後手動加上 `-lstdc++`。更常見的 C++ convention 是直接用：

```makefile
g++
```

或：

```makefile
c++
```

## 重要注意事項與可改進處

### 1. `copymatrix` 有維度 bug

目前：

```cpp
for(j=0;j<n;j++)
```

應該改成：

```cpp
for(j=0;j<p;j++)
```

否則非方陣會錯。

### 2. `matrices.h` 建議加入 include guard

可避免重複 include。

### 3. `readmatrix` 建議檢查 `fscanf` 回傳值

目前假設檔案格式一定正確。如果檔案數字不夠或格式錯，程式可能繼續用未預期的值。

### 4. `logdet` 對正定矩陣可用 Cholesky 更穩定

目前用 eigenvalue decomposition，可以工作，但對稱正定矩陣通常用 Cholesky 更有效也更穩定。

### 5. `inverse` 會修改原矩陣

這是刻意設計，但使用者要知道。如果還需要原矩陣，必須先 copy。

### 6. C 風格和 C++ 風格混用

這份程式混用了：

- C++ memory allocation：`new[]` / `delete[]`
- C memory allocation：`malloc` / `free`
- C file I/O：`FILE*`

在小型數值作業中很常見，但大型 C++ 專案會更偏好 `std::vector<double>`、RAII、`std::ifstream` / `std::ofstream`，減少手動管理記憶體的風險。

## 快速記憶版

如果只想抓重點：

- `matrices.h`：放函式宣告，讓其他檔案知道有哪些矩陣工具可以用。
- `matrices.cpp`：放函式實作，真正配置矩陣、讀寫矩陣、做矩陣運算。
- `double**`：這份程式用來表示二維矩陣。
- `allocmatrix` 配置後一定要用 `freematrix` 釋放。
- `matrixproduct` 是一般矩陣乘法。
- `dotmatrixproduct` 是逐元素相乘。
- `transposematrix` 會配置新矩陣並回傳，呼叫者要釋放。
- `inverse` 使用 LAPACKE 解 `A X = I` 來得到 `A^{-1}`，而且會覆蓋原矩陣。
- `logdet` 用 eigenvalues 估計 `log(det(A))`。
- 程式裡最明顯的 bug 是 `copymatrix` 內層應該用 `j < p`，不是 `j < n`。
