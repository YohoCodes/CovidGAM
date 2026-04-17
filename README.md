# CovidGAM: Smooth Deconvolution of Covid-19 Deaths

**CovidGAM** is an R-based statistical tool designed to infer the daily trajectory of new SARS-CoV-2 infections by deconvolving daily hospital death counts. The model uses a Penalized Poisson framework to account for the stochastic nature of death counts and the known delay between infection and death.

## Statistical Model

The script fits a **Penalized Poisson Deconvolution** model where daily deaths ($y_i$) on day $t_i$ are modeled as:

$$y_i \sim \text{Poi}(\mu_i)$$

The expected deaths $\mu_i$ are calculated as a convolution of the infection curve $f(t)$ and the probability function for the interval from infection to death $\pi(j)$:

$$\mu_i = \sum_{j=1}^{\min(29+i, 80)} f(t_i - j)\pi(j)$$

### Key Components:
* **Infection Curve ($f$):** Represented using **B-splines** ($f(t) = \sum b_k(t)\beta_k$) to allow for a flexible, smooth functional form.
* **Positivity Constraint:** Parameters are re-expressed as $\beta_k = \exp(\gamma_k)$ to ensure estimated infections remain strictly positive.
* **Smoothing Penalty:** To prevent overfitting, the model minimizes a penalized negative log-likelihood using a squared difference penalty on the spline coefficients:

$$P = \frac{\lambda}{2} \sum_{k=2}^{K-1}(\beta_{k-1} - 2\beta_k + \beta_{k+1})^2$$

* **Parameter Selection:** The optimal smoothing parameter $\lambda$ is chosen by minimizing the **Bayesian Information Criterion (BIC)**.


---

## Requirements

* **R** (Base R functions `optim`, `splines`, and `stats` are used for the core engine)
* **ggplot2** (Used for final visualization)

```r
install.packages("ggplot2")
```

---

## Files and Data

The script expects a data file (default: `engcov.txt`) in the working directory. The data must include:
* `julian`: Day of the year 2020.
* `nhs`: Daily death counts in English hospitals.

---

## How to Run

1.  **Set Working Directory:** Ensure your R session is set to the folder containing `CovidGAM.r` and `engcov.txt`.
2.  **Source the Script:**  
    * In **RStudio**, click **Source** or press `Ctrl+Shift+S`.  
    * In **VS Code**, run the file via the R extension's "Run Source" command.
3.  **Optimization:** The script will perform a grid search over $\log \lambda$ values (from $10^{-13}$ to $10^{-7}$) to find the BIC minimum.

---

## Interpreting Outputs

The script generates a final plot showing:
1.  **Actual vs. Fitted Deaths:** Daily death data overlaid with the model’s smoothed fit.
2.  **Estimated Infection Rate ($f$):** The deconvolved daily infection curve.
3.  **95% Confidence Intervals:** Generated via **non-parametric bootstrapping** (200 replicates) to assess the uncertainty of the infection trajectory.

## Performance Note
The full execution, including the BIC grid search and 200 bootstrap refits, is optimized for efficiency and should typically complete within a few dozen seconds. If it takes significantly longer, consider reducing the bootstrap count (`nb`) for testing.