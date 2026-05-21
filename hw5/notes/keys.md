# HW5 作業撰寫步驟

這份作業要完成兩個版本的 marginal likelihood 計算。Problem 1 使用一般 `double**` matrix helper functions；Problem 2 使用 GSL。兩題都要算同一個 regression `[1|A]`，測試資料是 `erdata.txt`，其中 `n = 158`、`p = 51`，第 1 欄是 response `D1`，其餘欄是 explanatory variables。

最後測試範例都應該印出：

```text
Marginal likelihood of regression [1|2,5,10] = -59.97893
```

注意：雖然 PDF 公式寫的是 marginal likelihood 本身，但範例輸出是負數，所以程式應該回傳 **log marginal likelihood**。

## 一、先理解要算什麼

給定變數集合 `A`，令 `k = lenA = |A|`。

你要建立：

- `D1`: data 第 1 欄，大小是 `n x 1`
- `DA`: data 中 A 指定的 explanatory columns，大小是 `n x k`
- `M = I_k + DA^T DA`，大小是 `k x k`

PDF 公式 (3) 的 log 版本可以寫成：

```text
log_marglik =
  lgamma((n + k + 2) / 2)
  - lgamma((k + 2) / 2)
  - 0.5 * log(det(M))
  - ((n + k + 2) / 2) * log(1 + D1^T D1 - D1^T DA M^{-1} DA^T D1)
```

在 C/C++ 裡要用 `lgamma()`、`log()`。

很重要的 index 規則：

- PDF 和 `A[] = {2,5,10}` 用的是 1-based column numbers。
- C/C++ array 是 0-based。
- 所以 response `D1` 是 `data[i][0]`。
- `A[j]` 對應到 `data[i][A[j] - 1]`。

## 二、Problem 1: `Matrices` 版本要做什麼

你要寫：

```cpp
double marglik(int n, int p, double** data, int lenA, int* A);
```

建議做法：

1. 在 `matrices.h` 加上 `marglik` function declaration。
2. 在 `matrices.cpp` 實作 `marglik`。
3. 建立一個 `main.cpp`，內容照 PDF 給的 main function。
4. 確認 makefile 會編譯 `main.cpp` 和 `matrices.cpp`，產生 executable `matrices`。

`marglik` 裡面的矩陣步驟：

1. 配置 `D1 = allocmatrix(n, 1)`。
2. 配置 `DA = allocmatrix(n, lenA)`。
3. 把 `data[i][0]` 複製到 `D1[i][0]`。
4. 把 `data[i][A[j] - 1]` 複製到 `DA[i][j]`。
5. 算 `DA^T`。
6. 算 `DA^T DA`。
7. 建立 `M = I + DA^T DA`。
8. 算 `logdet(M)`。要在 inverse 之前算，因為 inverse 會改掉矩陣。
9. 複製 `M` 到另一個矩陣 `Minv`，對 `Minv` 呼叫 `inverse(lenA, Minv)`。
10. 算 `D1^T D1`，這是一個 scalar。
11. 算 `DA^T D1`，大小是 `lenA x 1`。
12. 算 quadratic term: `(DA^T D1)^T Minv (DA^T D1)`。
13. 套入 log marginal likelihood 公式。
14. 釋放所有用 `allocmatrix` 配置的矩陣。

可以避免建立太多小矩陣的地方：

- `D1^T D1` 可以直接用 loop 加總 `data[i][0] * data[i][0]`。
- `DA^T D1` 可以用 loop 算成 `lenA x 1` vector。
- quadratic term 可以用雙層 loop：

```cpp
quad = 0;
for (i = 0; i < lenA; i++) {
  for (j = 0; j < lenA; j++) {
    quad += v[i] * Minv[i][j] * v[j];
  }
}
```

其中 `v[i] = DA^T D1`。

## 三、Problem 2: `MatricesGSL` 版本要做什麼

你要寫：

```cpp
double marglik(gsl_matrix* data, int lenA, int* A);
```

這題要求 relevant numerical functions 要來自 GSL，所以不要自己用 Problem 1 的 `allocmatrix`、`matrixproduct`、`inverse` 來做主要線性代數。

建議做法：

1. 建立 `MatricesGSL` 資料夾。
2. 放入 GSL 版本的 `matrices.cpp`、`matrices.h`、`main.cpp`、`makefile`。
3. `main.cpp` 使用 PDF 給的 GSL main function。
4. 在 `matrices.h` include GSL headers，例如：

```cpp
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
```

5. 在 `matrices.cpp` 實作 GSL 版本的 `marglik`。

GSL 版本建議矩陣步驟：

1. `n = data->size1`，`p = data->size2`，`k = lenA`。
2. 配置 `DA = gsl_matrix_alloc(n, k)`。
3. 配置 `M = gsl_matrix_alloc(k, k)`。
4. 配置 vector `d1 = gsl_vector_alloc(n)`。
5. 把 `gsl_matrix_get(data, i, 0)` 放進 `d1`。
6. 把 `gsl_matrix_get(data, i, A[j] - 1)` 放進 `DA(i,j)`。
7. 用 `gsl_blas_ddot(d1, d1, &yty)` 算 `D1^T D1`。
8. 用 `gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, DA, DA, 0.0, M)` 算 `DA^T DA`。
9. 對 `M` 加 identity。
10. 對 `M` 做 Cholesky decomposition，因為 `M = I + DA^T DA` 是 symmetric positive definite。
11. `logdet(M) = 2 * sum(log(diagonal entries of Cholesky factor))`。
12. 用 `gsl_blas_dgemv(CblasTrans, 1.0, DA, d1, 0.0, v)` 算 `v = DA^T D1`。
13. 解 `M x = v`，建議用 `gsl_linalg_cholesky_solve()`。
14. 用 `gsl_blas_ddot(v, x, &quad)` 算 `v^T M^{-1} v`。
15. 套入 log marginal likelihood 公式。
16. 釋放所有 GSL matrix/vector。

這個做法比直接算 inverse 更穩定，也比較符合 GSL 的使用方式。

## 四、資料夾與提交格式

根據 `requirements.md`，最後 zip 檔案要有兩個資料夾：

```text
Matrices/
  main.cpp
  matrices.cpp
  matrices.h
  makefile
  erdata.txt

MatricesGSL/
  main.cpp
  matrices.cpp
  matrices.h
  makefile
  erdata.txt
```

老師說提供的 extra `.txt` files 不要放，但你的程式執行需要 `erdata.txt`。如果課程要求資料檔另外提供，照課程規定；如果助教會直接在資料夾內執行，建議保留 `erdata.txt`，否則 `./matrices` 會找不到資料。

每個資料夾的 makefile 都要支援：

```bash
make clean
make all
```

而且 `make all` 後要產生 executable：

```text
matrices
```

## 五、在 cluster 上檢查

提交前一定要在 cluster 測：

```bash
module load GSL
module load R
make clean
make all
./matrices
```

Problem 1 如果使用 OpenBLAS/LAPACKE，也要確認 cluster 上有對應 module，例如：

```bash
module load OpenBLAS
```

如果 makefile 用了 `${EBROOTOPENBLAS}` 或 `${EBROOTGSL}`，要先確認 module load 後這些環境變數真的存在。

## 六、常見錯誤檢查

- 忘記 `A[j] - 1`，導致拿錯欄位。
- 把 response column 寫成 `data[i][1]`，正確是 `data[i][0]`。
- `copymatrix` 裡內層 loop 應該跑 `p`，不是 `n`。
- `logdet(M)` 要在 `M` 被 inverse 或 Cholesky overwrite 前算，或先 copy 一份。
- 忘記 free memory。
- makefile 沒有 `all` 或 `clean`。
- executable 名字不是 `matrices`。
- Problem 2 混用自己寫的 matrix inverse，而不是 GSL routine。
- 本機可以 compile，但 cluster 上缺 module 或 library path。

## 七、建議完成順序

1. 先完成 Problem 1 的 `main.cpp`、`marglik` declaration、`marglik` implementation。
2. 用 `make clean && make all && ./matrices` 確認輸出 `-59.97893`。
3. 修正 Problem 1 的 memory leak 和 index 問題。
4. 複製出 `MatricesGSL` 架構。
5. 改成 GSL data structure 和 GSL linear algebra。
6. 再次確認 GSL 版本輸出同樣是 `-59.97893`。
7. 整理資料夾，只保留需要提交的檔案。
8. 到 cluster 上做最後 compile/run 測試。
