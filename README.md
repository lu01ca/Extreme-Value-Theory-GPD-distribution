# Extreme Value Theory: GPD Modeling for Non-Stationary & Dependent Data

## Overview
This project explores **Extreme Value Theory (EVT)**, specifically the **Peaks Over Threshold (POT)** approach, to model the tail behavior of minimum daily temperature in Milan (1973-2024). 

Unlike simplified models that assume independent and identically distributed (i.i.d.) data, this analysis focuses on the rigorous statistical validation required for **non-stationary** and **serially dependent** environmental data, ensuring that the theoretical assumptions for the Generalized Pareto Distribution (GPD) are strictly met.

## Theoretical Background

### From GEV to GPD
The foundation of EVT relies on the **Generalized Extreme Value (GEV)** distribution, which models block maxima and is defined as:
$$G(z) = \exp\left\{-\left[1 + \xi \left(\frac{z-\mu}{\sigma}\right)\right]^{-1/\xi}\right\}$$
However, the GEV approach wastes data by only considering maximums over large blocks. The **POT approach** improves efficiency by modeling all exceedances above a sufficiently high threshold $u$. Asymptotically, these exceedances follow a **Generalized Pareto Distribution (GPD)**:
$$H(y) = 1 - \left(1 + \frac{\xi y}{\tilde{\sigma}}\right)^{-1/\xi}$$
where $y = (x - u) > 0$ are the threshold excesses, and $\tilde{\sigma} = \sigma + \xi(u-\mu)$ is the scale parameter.

### The Shape Parameter ($\xi$)
The shape parameter $\xi$ is crucial as it dictates the behavior of the distribution's tail:
* **$\xi > 0$:** Heavy, unbounded tail (indicates a high probability of extreme events).
* **$\xi < 0$:** Bounded upper tail (implies a finite maximum possible limit).
* **$\xi = 0$:** Light, exponential tail.

## Statistical Pre-processing
Environmental data like daily minimum temperature rarely satisfy the i.i.d. assumption. Before fitting the GPD, the data underwent a comprehensive diagnostic and transformation phase:

### 1. Handling Non-Stationarity (Seasonality)
Initial time-series plotting and formal testing using the **Augmented Dickey-Fuller (ADF) Test** revealed significant seasonal non-stationarity. To obtain identically distributed subsets, the dataset was split into **four meteorological seasons** (Winter, Spring, Summer, Autumn). Each season was modeled independently to capture specific atmospheric tail behaviors.

### 2. Handling Dependency (Declustering)
temperature exhibit strong serial correlation (e.g., a temperature day is often followed by another). Since EVT requires independent exceedances, a **Declustering** technique using the *Runs Method* was applied. 
This method filters the data by grouping consecutive exceedances into clusters, separated by a run length of $r$ consecutive observations below the threshold $u$. Only the maximum value of each cluster is extracted, ensuring the independence assumption for the GPD is valid.

## Modeling Methodology
The core of the analysis involves the following steps implemented in R for each season:

1.  **Threshold Selection ($u$)**: Utilization of the **Mean Residual Life Plot** (`mrlplot`) to identify the optimal threshold where the mean excess becomes linear.
2.  **Parameter Estimation**: Fitting the GPD using **Maximum Likelihood Estimation (MLE)** to determine the shape ($\xi$) and scale ($\sigma$) parameters.
3.  **Model Diagnostics**: Verification of the goodness-of-fit through **Probability Plots**, **Quantile Plots**, **Density Plots**, and **Return Level Plots** to estimate $N$-year return periods.

## Application & Results
*(Note: Insert here a brief summary of your results. Example: "The analysis revealed that Winter exhibits a heavier tail ($\xi > 0$) compared to Summer ($\xi < 0$), indicating a higher risk of extreme wind gusts during the colder months.")*

![Diagnostic Plots](link-to-your-image.png)
*(Note: Add an image of your `gpd.diag` output or `mrlplot` here)*

## Repository Structure
* `Latex.R`: Complete R script including EDA, stationarity/dependency tests, declustering, and model fitting.
* `EVT_Pareto_generalizzata.pdf`: Detailed technical report (in Italian) covering the mathematical framework, proofs, and numerical results.

## Requirements & Environment
The analysis was performed in R and requires the following libraries:
```R
library(ismev)   # Extreme value modeling
library(evd)     # Extreme value distributions functions
library(tseries) # For ADF Test
```
## References & License
Methodology Reference: Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values. Springer.
