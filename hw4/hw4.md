# HW4 程式說明

這份作業要做的是「一個一個變數」的 logistic regression。資料有 61 欄：前 60 欄是解釋變數 `X1` 到 `X60`，最後一欄是二元反應變數 `Y`，也就是每筆資料的結果只有 `0` 或 `1`。

程式檔案是 `hw4-solution.R`。它完成三件事：

1. 用 Newton-Raphson 找 Bayesian logistic regression 的 posterior mode。
2. 用 Laplace approximation 算每個模型的 log marginal likelihood。
3. 用 Metropolis-Hastings 抽樣，估計 `beta0` 和 `beta1` 的 posterior mean。

## Logistic Regression 在做什麼

對每一個解釋變數 `x`，我們建立模型：

```text
log(P(Y = 1 | x) / P(Y = 0 | x)) = beta0 + beta1 * x
```

`beta0` 是截距，`beta1` 是這個變數的效果。如果 `beta1` 很大，表示 `x` 增加時，`Y = 1` 的機率會明顯改變。

程式中用 `inv_logit()` 把線性值 `eta = beta0 + beta1*x` 轉成機率：

```text
pi = exp(eta) / (1 + exp(eta))
```

也就是 `pi = P(Y = 1 | x)`。

## Prior 和 Posterior

作業假設：

```text
beta0 ~ N(0, 1)
beta1 ~ N(0, 1)
```

這代表在看資料以前，我們相信 `beta0` 和 `beta1` 大概在 0 附近，但資料仍然可以把它們推到其他值。

程式的 `log_posterior_kernel()` 算的是作業 PDF 裡的 `l*(beta0, beta1)`：

```text
l* = log prior + log likelihood
```

我們用 log scale 計算，因為很多機率相乘會非常小，直接乘容易在電腦裡變成 0。

## Newton-Raphson

`newton_raphson_mode()` 用來找 posterior distribution 最高的地方，也就是 posterior mode：

```text
beta_new = beta_old - Hessian^{-1} * Gradient
```

在程式裡：

```r
grad <- gradient_log_posterior(beta, x, y)
hess <- hessian_log_posterior(beta, x, y)
step <- solve(hess, grad)
candidate <- beta - step
```

直覺上，gradient 告訴我們「往哪個方向會增加」，Hessian 告訴我們「曲線彎曲程度」。Newton-Raphson 用這兩個資訊很快地跳到最高點附近。

程式也有一個小保護：如果某次 Newton step 太大，造成數值不穩或讓 `l*` 變差，就把步伐縮小。這樣對不同 predictor 比較穩。

## Laplace Approximation

Problem 1 要算 marginal likelihood：

```text
P(D) = integral exp(l*(beta0, beta1)) d beta0 d beta1
```

這個積分通常不能直接算，所以用 Laplace approximation。核心想法是：posterior 在 mode 附近最重要，所以用一個二維常態分布去近似 mode 附近的形狀。

程式的 `laplace_log_marginal_likelihood()` 算的是：

```text
log P(D) approx =
log(2*pi) + l*(mode) - 0.5 * log(det(-Hessian at mode))
```

注意我們回傳的是 log marginal likelihood，不是 marginal likelihood 本身。這也是 PDF 特別要求的，因為原本的值可能非常小。

## Metropolis-Hastings

Problem 2 要從 posterior distribution 抽 10000 個 `(beta0, beta1)` 樣本，然後取平均：

```text
beta0_bayes = average of sampled beta0
beta1_bayes = average of sampled beta1
```

程式的 `metropolis_hastings()` 做這件事。

流程是：

1. 從 posterior mode 開始。
2. 用二維常態分布提出一個新的候選點。
3. 如果候選點的 `l*` 比目前點更大，就接受。
4. 如果候選點的 `l*` 比目前點更小，也不是完全拒絕，而是用機率接受。

接受規則用 log scale 寫成：

```r
log(runif(1)) <= proposed_log_post - current_log_post
```

這和 PDF 裡的接受機率：

```text
min(1, exp(l*(proposal) - l*(current)))
```

是同一件事，只是 log scale 比較穩。

## `bayesLogistic()` 做什麼

`bayesLogistic(apredictor, response, data, NumberOfIterations)` 是作業 starter code 要我們完成的主函數。

它會針對某一個 predictor 做完整分析：

1. 取出第 `apredictor` 欄當作 `x`。
2. 取出最後一欄當作 `y`。
3. 用 Newton-Raphson 找 posterior mode。
4. 用 Laplace approximation 算 log marginal likelihood。
5. 用 Metropolis-Hastings 算 Bayesian estimates。
6. 用 `glm()` 算 MLE，當作檢查用的比較值。

最後回傳一個 list，包含：

```text
apredictor
logmarglik
beta0bayes
beta1bayes
beta0mle
beta1mle
beta0mode
beta1mode
acceptance_rate
```

## `main()` 做什麼

`main()` 是 Problem 3 的平行版本。

它讀入資料後，會對 60 個 predictor 分別呼叫 `bayesLogistic()`：

```r
results <- clusterApply(cluster, 1:lastPredictor, bayesLogistic,
                        response, data, NumberOfIterations)
```

如果電腦有安裝 `snow`，程式會用 `snow`。如果沒有安裝，程式會自動改用 R 內建的 `parallel` package，因為這台機器目前沒有 `snow`。

每一列輸出代表一個 regression，例如：

```text
Regression of Y on explanatory variable 3 has log marginal likelihood ...
```

括號外的 `beta0` 和 `beta1` 是 Bayesian posterior mean，括號內的是 `glm()` 算出的 MLE。

## 如何執行

在 terminal 裡：

```bash
Rscript hw4-solution.R
```

它會自動找資料檔。PDF 寫的是 `534binarydata.txt`，但這個資料夾裡的檔名是 `534binarydata-1.txt`，所以程式有處理這個差異。

如果你只想快速測試一小段，可以進 R 後執行：

```r
source("hw4-solution.R")
results <- main("534binarydata-1.txt", NumberOfIterations = 100, clusterSize = 1)
```

正式作業使用 10000 iterations：

```r
results <- main("534binarydata-1.txt", NumberOfIterations = 10000, clusterSize = 10)
```

## 檔案中重要函數總表

`log1pexp()`：穩定計算 `log(1 + exp(x))`。

`inv_logit()`：把任何實數轉成 0 到 1 之間的機率。

`log_likelihood()`：計算 logistic regression 的 log-likelihood。

`log_posterior_kernel()`：計算 PDF 裡的 `l*`。

`gradient_log_posterior()`：計算 `l*` 的 gradient。

`hessian_log_posterior()`：計算 `l*` 的 Hessian。

`newton_raphson_mode()`：用 Newton-Raphson 找 posterior mode。

`laplace_log_marginal_likelihood()`：用 Laplace approximation 算 log marginal likelihood。

`metropolis_hastings()`：用 Metropolis-Hastings 抽 posterior samples 並估計 beta。

`bayesLogistic()`：針對單一 predictor 跑完整 Bayesian logistic regression。

`main()`：針對全部 60 個 predictors 跑完整作業。
