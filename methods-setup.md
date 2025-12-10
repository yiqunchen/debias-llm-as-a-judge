# PPI vs Rogan–Gladen (RG) in Binary Misclassification

We compare two estimators of a scalar parameter
\(\theta = \mathbb{E}[Z] = \Pr(Z=1)\):

1. Rogan–Gladen (RG) misclassification‐corrected estimator  
2. Vanilla Prediction–Powered Inference (PPI)

and then show how both fit into a more general “PPI with transform \(g\)” framework.

---

## 1. Setup and notation (binary misclassification)

- True (human) label: \(Z \in \{0,1\}\)
- Target parameter: \(\theta = \mathbb{E}[Z] = \Pr(Z=1)\)
- Surrogate / judge label: \(\hat Z \in \{0,1\}\)

Define:

- Sensitivity:
  \[
  q_1 = \Pr(\hat Z = 1 \mid Z = 1)
  \]
- Specificity:
  \[
  q_0 = \Pr(\hat Z = 0 \mid Z = 0)
  \]
- Youden’s \(J\):
  \[
  J = q_0 + q_1 - 1 \in (0,1] \quad\text{(assume the judge is better than random)}
  \]

### Data

- **Big test set** of size \(n\): only \(\hat Z_i\), \(i=1,\dots,n\)
  \[
  \hat p = \frac{1}{n} \sum_{i=1}^n \hat Z_i \approx p := \Pr(\hat Z = 1).
  \]

- **Small calibration set** of size \(m\): both \(Z_j,\hat Z_j\), \(j=1,\dots,m\).

**Goal**: estimate \(\theta = \mathbb{E}[Z]\).

---

## 2. Two estimators

### 2.1 Rogan–Gladen (RG) estimator

Under the standard misclassification model,
\[
p = \Pr(\hat Z = 1)
  = q_1 \theta + (1-q_0)(1-\theta)
  = J\,\theta + (1-q_0).
\]

Solve for \(\theta\):
\[
\theta = \frac{p + q_0 - 1}{J}.
\]

Plug in empirical estimates:
- \(\hat p = \frac{1}{n}\sum_{i=1}^n \hat Z_i\)
- \(\hat q_0 = \frac{\#\{j: \hat Z_j=0, Z_j=0\}}{\#\{j: Z_j=0\}}\)
- \(\hat q_1 = \frac{\#\{j: \hat Z_j=1, Z_j=1\}}{\#\{j: Z_j=1\}}\)

The RG estimator is:
\[
\hat\theta_{\text{RG}}
  = \frac{\hat p + \hat q_0 - 1}{\hat q_0 + \hat q_1 - 1}.
\]

If we pretend \(q_0,q_1\) are known (or \(m\) is huge so their uncertainty is negligible), this is just a linear transform of \(\hat p\):
\[
\hat\theta_{\text{RG}} \approx \frac{\hat p + q_0 - 1}{J}.
\]

---

### 2.2 Vanilla PPI estimator

Let:
- \(Y := Z\)
- \(f(X) := \hat Z\)  (we ignore \(X\) explicitly; treat \(\hat Z\) as a function of data)

Then:
\[
\theta
  = \mathbb{E}[Z]
  = \mathbb{E}[\hat Z] - \mathbb{E}[\hat Z - Z]
  = p - \Delta,
\]
where
\[
\Delta := \mathbb{E}[\hat Z - Z].
\]

Estimators:

- On the test set:
  \[
  \hat p = \frac{1}{n} \sum_{i=1}^n \hat Z_i.
  \]

- On the calibration set:
  \[
  \hat\Delta
    = \frac{1}{m} \sum_{j=1}^m (\hat Z_j - Z_j).
  \]

PPI estimator:
\[
\hat\theta_{\text{PPI}} = \hat p - \hat\Delta.
\]

This uses only the *average prediction error* from the calibration set; it never explicitly uses \(q_0,q_1\).

---

## 3. Asymptotic variances (simple regime)

Assumptions for this comparison:

- Test and calibration sets are i.i.d. from the same population.
- \(n\) (test size) and \(m\) (calibration size) are both large.
- Treat \((q_0,q_1)\), and hence \(J\), as fixed constants (ignoring their estimation variance for now).

### 3.1 Variance of RG estimator (known \(q_0,q_1\))

If \(q_0,q_1\) are known, then
\[
\hat\theta_{\text{RG}} = \frac{\hat p + q_0 - 1}{J}
\]
is linear in \(\hat p\).

Since \(\hat p\) is a binomial proportion:
\[
\mathrm{Var}(\hat p) \approx \frac{p(1-p)}{n},
\]
we obtain:
\[
\mathrm{Var}(\hat\theta_{\text{RG}})
  \approx \frac{p(1-p)}{n J^2}.
\tag{RG}
\]

This is the usual binomial variance scaled by \(1/J^2\).

---

### 3.2 Variance of PPI estimator

General result for mean PPI (Angelopoulos et al.):

If \(\hat\theta\) is constructed from:
- a large “unlabeled” set of size \(n_{\text{unlab}}\) with only \(f(X)\),
- a “labeled” set of size \(n_{\text{lab}}\) with both \(f(X)\) and \(Y\),

then asymptotically:
\[
\mathrm{Var}(\hat\theta_{\text{PPI}})
  \approx
  \frac{\mathrm{Var}(f(X))}{n_{\text{unlab}}}
  +
  \frac{\mathrm{Var}(f(X) - Y)}{n_{\text{lab}}}.
\]

In our case:
- \(f(X) = \hat Z\),
- \(Y = Z\),
- \(n_{\text{unlab}}=n\) (test set),
- \(n_{\text{lab}}=m\) (calibration set),

so:
\[
\mathrm{Var}(\hat\theta_{\text{PPI}})
  \approx
  \frac{\mathrm{Var}(\hat Z)}{n}
  +
  \frac{\mathrm{Var}(\hat Z - Z)}{m}.
\tag{PPI}
\]

We have:
\[
\mathrm{Var}(\hat Z) = p(1-p).
\]

Also, \(\hat Z - Z \in \{-1,0,1\}\):

- \(\hat Z - Z = +1\) for false positives (FP),
- \(\hat Z - Z = -1\) for false negatives (FN),
- \(\hat Z - Z = 0\) otherwise.

Let:
- \(\Pr(\text{FP}) = \Pr(\hat Z = 1, Z = 0)\),
- \(\Pr(\text{FN}) = \Pr(\hat Z = 0, Z = 1)\).

Let \(\Delta = \mathbb{E}[\hat Z - Z] = p - \theta\). Then:
\[
\mathrm{Var}(\hat Z - Z)
  = \Pr(\text{FP}) + \Pr(\text{FN}) - \Delta^2.
\]

Therefore:
\[
\mathrm{Var}(\hat\theta_{\text{PPI}})
  \approx
  \frac{p(1-p)}{n}
  +
  \frac{\Pr(\text{FP}) + \Pr(\text{FN}) - \Delta^2}{m}.
\]

If \(m\) is large enough that the calibration term is negligible:
\[
\mathrm{Var}(\hat\theta_{\text{PPI}})
  \approx \frac{p(1-p)}{n}.
\]

---

### 3.3 Direct comparison when calibration is large

In the “large calibration” regime (ignore the calibration term in (PPI)):

- PPI variance:
  \[
  \mathrm{Var}(\hat\theta_{\text{PPI}})
    \approx \frac{p(1-p)}{n}.
  \]

- RG variance:
  \[
  \mathrm{Var}(\hat\theta_{\text{RG}})
    \approx \frac{p(1-p)}{n J^2}.
  \]

Since \(0 < J \le 1\), we have \(1/J^2 \ge 1\), with strict inequality whenever the judge is imperfect (i.e., \(J < 1\)). Thus:
\[
\mathrm{Var}(\hat\theta_{\text{PPI}})
  \le
  \mathrm{Var}(\hat\theta_{\text{RG}}),
\]
with equality only if the judge is perfect (\(J = 1\)).

With **finite** \(m\), PPI gains the additional term
\(\mathrm{Var}(\hat Z - Z)/m\), and RG gains extra variance from estimating \(\hat q_0,\hat q_1\) (not written out here), so there is no uniform dominance in all parameter regimes. But the structural difference remains:

- RG amplifies test‑set noise by \(1/J^2\),
- PPI keeps the test‑set noise at the base \(p(1-p)/n\) and adds a separate calibration term.

---

## 4. PPI with a general transform \(g\)

Let \(S\) denote the surrogate (here \(S = \hat Z\)), and let \(X\) be any additional features (possibly ignored).

For **any** function \(g(S,X)\), define the PPI‑style estimator:
\[
\hat\theta_g
  =
  \underbrace{
    \frac{1}{n} \sum_{\text{test}} g(S_i, X_i)
  }_{\text{plug-in on big sample}}
  -
  \underbrace{
    \frac{1}{m} \sum_{\text{cal}} \big( g(S_j, X_j) - Z_j \big)
  }_{\text{rectifier on calibration}}.
\]

Then:
\[
\mathbb{E}[\hat\theta_g]
  = \mathbb{E}[g(S,X)] - \mathbb{E}[g(S,X) - Z]
  = \mathbb{E}[Z] = \theta.
\]

Asymptotic variance:
\[
\mathrm{Var}(\hat\theta_g)
  \approx
  \frac{\mathrm{Var}(g(S,X))}{n}
  +
  \frac{\mathrm{Var}(g(S,X) - Z)}{m}.
  \tag{*}
\]

Choosing \(g\) is, therefore, a variance‑optimization problem.

### 4.1 Special choices of \(g\)

#### (a) Identity / vanilla PPI

\[
g_{\text{id}}(S) = S = \hat Z.
\]

Then:
\[
\mathrm{Var}(\hat\theta_{\text{id}})
  \approx
  \frac{\mathrm{Var}(\hat Z)}{n}
  +
  \frac{\mathrm{Var}(\hat Z - Z)}{m}.
\]

#### (b) Rogan–Gladen transform as \(g\)

Suppose you have fixed or pilot estimates \(\tilde q_0,\tilde q_1\) of \((q_0,q_1)\), and define \(\tilde J = \tilde q_0 + \tilde q_1 - 1\).

Define:
\[
g_{\text{RG}}(S) = \frac{S + \tilde q_0 - 1}{\tilde J}.
\]

Then:
\[
\frac{1}{n} \sum_{\text{test}} g_{\text{RG}}(S_i)
  =
  \frac{\hat p + \tilde q_0 - 1}{\tilde J},
\]
which is the RG plug‑in estimator with \(\tilde q_0,\tilde q_1\) in place of \(\hat q_0,\hat q_1\).

If you **stop here** (omit the rectifier), you recover the usual RG‑style estimator (with fixed \(q_0,q_1\)).

With full PPI:
\[
\hat\theta_{g_{\text{RG}}}
  =
  \frac{1}{n} \sum_{\text{test}} g_{\text{RG}}(S_i)
  -
  \frac{1}{m} \sum_{\text{cal}}
    \big( g_{\text{RG}}(S_j) - Z_j \big),
\]
which is still unbiased for \(\theta\) and uses the calibration set to correct any residual bias in \(g_{\text{RG}}\).

#### (c) Learned / optimal \(g\)

More generally, you can **learn** \(g\) on the calibration set (e.g., via regression) so that \(g(S,X)\) approximates \(Z\) as well as possible, thereby minimizing the variance in \((*)\). RePPI shows how to construct such an asymptotically optimal \(g\) via influence‑function regression.

---

## 5. Summary

- **RG estimator**:
  \[
  \hat\theta_{\text{RG}}
    =
    \frac{\hat p + \hat q_0 - 1}{\hat q_0 + \hat q_1 - 1}.
  \]

- **Vanilla PPI estimator (with \(g(S)=S\))**:
  \[
  \hat\theta_{\text{PPI}}
    =
    \hat p - \hat\Delta
    =
    \frac{1}{n} \sum_{\text{test}} \hat Z_i
    -
    \frac{1}{m} \sum_{\text{cal}} (\hat Z_j - Z_j).
  \]

- **General PPI with transform \(g\)**:
  \[
  \hat\theta_g =
    \frac{1}{n} \sum_{\text{test}} g(S_i,X_i)
    -
    \frac{1}{m} \sum_{\text{cal}} (g(S_j,X_j) - Z_j).
  \]

In the large‑calibration limit and when the judge is imperfect (\(J<1\)), PPI has asymptotic variance \(\approx p(1-p)/n\), whereas RG has variance \(\approx p(1-p)/(nJ^2)\); thus, PPI is asymptotically more efficient in that regime.

---

# Multiclass (K-Class) Misclassification: RG vs PPI

We now consider a **K-class** categorical outcome, observed through a noisy surrogate (e.g., an LLM labeler or classifier). We compare:

- A **multiclass Rogan–Gladen / confusion-matrix inversion** estimator, and  
- A **Prediction-Powered Inference (PPI)** estimator,

and then show how both fit into a general PPI-with-transform-\(g\) framework.

## 1. Setup and Notation (K-Class Categorical Misclassification)

Let the true label be \(Y \in \{1,\dots,K\}\) with target class probabilities \(\boldsymbol\pi = (\pi_1,\dots,\pi_K)^\top\). The surrogate \(S\in\{1,\dots,K\}\) obeys a misclassification matrix \(M\) with entries \(M_{ab} = \Pr(S=a \mid Y=b)\) and observed surrogate probabilities \(\mathbf p = (p_1,\dots,p_K)^\top\) satisfying \(\mathbf p = M\boldsymbol\pi\) (assume \(M\) invertible).

Data:

- **Large test set** of size \(n\): observe \(S_i\) only; \(\hat{\mathbf p}\) is the empirical class distribution of \(S\).
- **Calibration set** of size \(m\): observe \((Y_j,S_j)\); estimate \(M\) via \(\hat M_{ab} = \#\{j:S_j=a,Y_j=b\} / \#\{j:Y_j=b\}\).

Goal: estimate \(\boldsymbol\pi\).

## 2. Two Estimators

### 2.1 Multiclass Rogan–Gladen / Confusion-Matrix Inversion

Since \(\mathbf p = M\boldsymbol\pi\), the true prevalence vector is \(\boldsymbol\pi = M^{-1}\mathbf p\). Plugging in \(\hat M\) and \(\hat{\mathbf p}\) yields:
\[
  \hat{\boldsymbol\pi}_{\text{RG}} = \hat M^{-1} \hat{\mathbf p}.
\]
Componentwise \(\hat\pi_{\text{RG},k}\) is the \(k\)-th entry of \(\hat M^{-1}\hat{\mathbf p}\). This generalizes binary RG by inverting the estimated confusion matrix.

### 2.2 PPI Estimator (Identity \(g\), Class-by-Class)

Represent \(Y\) via one-hot vectors \(\mathbf Y = (1\{Y=1\},\dots,1\{Y=K\})^\top\), so \(\boldsymbol\pi = \mathbb{E}[\mathbf Y]\). For class \(k\), set \(g^{\text{id}}_k(S) = 1\{S=k\}\) and apply PPI:
\[
  \hat\pi^{\text{PPI}}_k = \frac{1}{n}\sum_{i=1}^n g^{\text{id}}_k(S_i) + \frac{1}{m}\sum_{j=1}^m \big(Y_j^{(k)} - g^{\text{id}}_k(S_j)\big).
\]
Stacking gives:
\[
  \hat{\boldsymbol\pi}_{\text{PPI}} = \frac{1}{n}\sum_{i=1}^n g^{\text{id}}(S_i) + \frac{1}{m}\sum_{j=1}^m (\mathbf Y_j - g^{\text{id}}(S_j)).
\]

## 3. Asymptotic Variances (Sketch)

Treat estimators as vectors in \(\mathbb{R}^K\).

### 3.1 RG Estimator (Known \(M\))

If \(M\) is known, \(\hat{\boldsymbol\pi}_{\text{RG}} = M^{-1}\hat{\mathbf p}\) with multinomial covariance:
\[
  \mathrm{Cov}(\hat{\boldsymbol\pi}_{\text{RG}}) \approx \frac{1}{n}\, M^{-1}\big(\mathrm{diag}(\mathbf p) - \mathbf p\mathbf p^\top\big)(M^{-1})^\top.
\]
Estimating \(M\) adds a delta-method variance term.

### 3.2 PPI Estimator (Vector Form)

For any vector-valued \(g(S)\):
\[
  \hat{\boldsymbol\pi}_g = \frac{1}{n}\sum_{i=1}^n g(S_i) + \frac{1}{m}\sum_{j=1}^m (\mathbf Y_j - g(S_j)),
\]
with
\[
  \mathrm{Cov}(\hat{\boldsymbol\pi}_g) \approx \frac{1}{n}\,\mathrm{Cov}(g(S)) + \frac{1}{m}\,\mathrm{Cov}(\mathbf Y - g(S)).
\]
Setting \(g = g^{\text{id}}\) produces the classwise PPI variances. Other \(g\) choices target lower variance.

## 4. PPI with Transform \(g\) (Multiclass)

Let \(g(S,X)\in\mathbb{R}^K\). Then
\[
  \hat{\boldsymbol\pi}_g = \frac{1}{n}\sum_{i=1}^n g(S_i,X_i) + \frac{1}{m}\sum_{j=1}^m (\mathbf Y_j - g(S_j,X_j))
\]
remains unbiased with covariance \(\frac{1}{n}\mathrm{Cov}(g) + \frac{1}{m}\mathrm{Cov}(\mathbf Y-g)\). Examples:

1. **Identity / naïve PPI**: \(g^{\text{id}}(S)\) (one-hot).  
2. **Confusion-matrix transform (RG)**: with \(W=\hat M^{-1}\), set \(g^{\text{RG}}(S)=W e_S\) (column \(S\) of \(W\)). The plug-in term alone recover RG; adding the rectifier keeps unbiasedness while reducing variance.  
3. **Soft/learned \(g\)**: fit \(g_k(S,X)\approx \Pr(Y=k\mid S,X)\), e.g., via multinomial logistic regression, and plug into the same formula.

This shows RG is a special linear \(g\), while PPI can incorporate RG or learned transforms to trade off bias and variance flexibly.

### 4.2 Projection onto the simplex

Rather than ad-hoc clipping of each coordinate, we can enforce the probability-simplex constraints on any multiclass estimator by projecting:
\[
  \tilde{\boldsymbol\pi}
    = \arg\min_{\boldsymbol\pi \in \Delta_{K-1}} \|\boldsymbol\pi - \hat{\boldsymbol\pi}\|_2^2,
\qquad
  \Delta_{K-1} = \{\boldsymbol\pi: \pi_k \ge 0,\ \sum_{k=1}^K \pi_k = 1\}.
\]
The minimizer has a closed-form “sorting + threshold” solution identical to the sparsemax/simplex-projection algorithm. It returns the closest valid probability vector to the unconstrained estimate and is standard in quantification and survey calibration workflows.

### 4.3 Constrained RG / constrained likelihood

On the RG side, instead of computing \(\hat M^{-1} \hat{\mathbf p}\) directly, we can solve the constrained least-squares problem
\[
  \min_{\boldsymbol\pi \in \Delta_{K-1}} \|\hat{\mathbf p} - \hat M \boldsymbol\pi\|_2^2,
\]
or the constrained multinomial likelihood. This enforces non-negativity and the sum-to-one constraint “inside” the estimation, and can be viewed as RG with a simplex constraint baked in at the parameter level. The projection-based RG and the constrained least-squares RG coincide when \(\hat M\) is invertible and the solution lies in the simplex interior, but differ when inversion would leave the simplex. Both are implemented in the companion code so you can compare unconstrained, projected, and fully constrained variants.
