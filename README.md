# Extreme Value Theory: GPD Modeling for Non-Stationary & Dependent Data

## Overview
This project explores **Extreme Value Theory (EVT)**, specifically the **Peaks Over Threshold (POT)** approach, to model the tail behavior of **daily minimum temperatures** in Milan (1973-2024). 

Unlike simplified models that assume independent and identically distributed (i.i.d.) data, this analysis focuses on the rigorous statistical validation required for **non-stationary** and **serially dependent** environmental data, ensuring that the theoretical assumptions for the Generalized Pareto Distribution (GPD) are strictly met.

*Note: Since standard EVT models maxima, the dataset was transformed using $X' = -X$ to analyze extreme minima (extreme cold events). Return levels were then re-transformed to their original scale.*

## Theoretical Background

### From GEV to GPD
The foundation of EVT relies on the **Generalized Extreme Value (GEV)** distribution, which models block maxima and is defined as:
$$G(z) = \exp\{ -[ 1 + \xi ( \frac{z-\mu}{\sigma} ) ]^{-1/\xi} \}$$
However, the GEV approach wastes data by only considering maximums over large blocks. The **POT approach** improves efficiency by modeling all exceedances above a sufficiently high threshold $u$. Asymptotically, these exceedances follow a **Generalized Pareto Distribution (GPD)**:
$$H(y) = 1 - \left(1 + \frac{\xi y}{\tilde{\sigma}}\right)^{-1/\xi}$$
where $y = (x - u) > 0$ are the threshold excesses, and $\tilde{\sigma} = \sigma + \xi(u-\mu)$ is the scale parameter.

### The Shape Parameter ($\xi$)
The shape parameter $\xi$ is crucial as it dictates the behavior of the distribution's tail:
* **$\xi > 0$:** Heavy, unbounded tail (indicates a high probability of extreme events).
* **$\xi < 0$:** Bounded upper tail (implies a finite maximum possible limit).
* **$\xi = 0$:** Light, exponential tail.

## Statistical Pre-processing
Environmental data like daily minimum temperatures rarely satisfy the i.i.d. assumption. Before fitting the GPD, the data underwent a comprehensive diagnostic and transformation phase:

### 1. Handling Non-Stationarity (Seasonality)
Initial time-series plotting and formal testing using the **Augmented Dickey-Fuller (ADF) Test** revealed significant seasonal non-stationarity. To obtain identically distributed subsets, the dataset was split into **four meteorological seasons** (Winter, Spring, Summer, Autumn). Each season was modeled independently to capture specific atmospheric tail behaviors.

### 2. Handling Dependency (Declustering)
Temperatures exhibit strong serial correlation (e.g., a cold day is often followed by another cold day). Since EVT requires independent exceedances, a **Declustering** technique using the *Runs Method* was applied. 
This method filters the data by grouping consecutive exceedances into clusters, separated by a run length of $r$ consecutive observations below the threshold $u$. Only the maximum value of each cluster is extracted, ensuring the independence assumption for the GPD is valid.

## Modeling Methodology
The core of the analysis involves the following steps implemented in R for each season:

1.  **Threshold Selection ($u$)**: Utilization of the **Mean Residual Life Plot** (`mrlplot`) to identify the optimal threshold where the mean excess becomes linear.
2.  **Parameter Estimation**: Fitting the GPD using **Maximum Likelihood Estimation (MLE)** to determine the shape ($\xi$) and scale ($\sigma$) parameters.
3.  **Model Diagnostics**: Verification of the goodness-of-fit through **Probability Plots**, **Quantile Plots**, **Density Plots**, and **Return Level Plots** to estimate $N$-year return periods.

![diagnostica.pdf](https://github.com/user-attachments/files/25427764/diagnostica.pdf)

## Application & Results: Return Levels for Extreme Cold
The model successfully estimated return levels ($\hat{z}_T$) for extreme cold events across different return periods ($T$), adjusted for temporal dependency ($\theta$) and computed using the Delta Method for 95% Confidence Intervals.

| Season | $T$ (years) | $\hat{z}_T$ (°C) | 95% CI |
| :--- | :--- | :--- | :--- |
| **DJF (Winter)** | 10 | -11.3 | [-12.6, -10.1] |
| | 20 | -12.8 | [-14.3, -11.2] |
| | 50 | -14.6 | [-16.5, -12.7] |
| | 100 | -16.0 | [-18.2, -13.8] |
| **MAM (Spring)** | 10 | -4.3 | [-5.4, -3.3] |
| | 20 | -5.4 | [-6.7, -4.1] |
| | 50 | -6.8 | [-8.4, -5.1] |
| | 100 | -7.8 | [-9.8, -5.9] |
| **JJA (Summer)** | 10 | +7.4 | [+6.3, +8.5] |
| | 20 | +6.2 | [+4.8, +7.5] |
| | 50 | +4.6 | [+2.8, +6.3] |
| | 100 | +3.3 | [+1.4, +5.3] |
| **SON (Autumn)** | 10 | -4.8 | [-5.4, -4.2] |
| | 20 | -5.4 | [-6.0, -4.7] |
| | 50 | -5.9 | [-6.7, -5.1] |
| | 100 | -6.2 | [-7.1, -5.2] |

### Key Insights
* **Winter (DJF):** Exhibits the most severe extremes. The 100-year return level of -16.0°C is highly consistent with the historical record of Milan Linate (-14.4°C in Jan 1985).
* **Spring (MAM):** While less intense than winter, a 100-year event of -7.8°C represents an exceptional thermal anomaly compared to the seasonal average.
* **Summer (JJA):** Return levels remain positive even for $T=100$ years (+3.3°C). While statistically correct, this has a limited physical interpretation: these are not "extreme cold" events in the traditional sense, but rather unusually cool nights for the summer period.
* **Autumn (SON):** The narrow gap between $T=10$ (-4.8°C) and $T=100$ (-6.2°C) is a direct consequence of a bounded upper tail ($\xi < 0$), reflecting that autumn temperatures, while dropping below zero, do not reach the extreme lows typical of winter.

## Repository Structure
* `code.R`: Complete R script including EDA, stationarity/dependency tests, declustering, and model fitting.
* `EVT_Pareto_generalizzata.pdf`: Detailed technical report (in Italian) covering the mathematical framework, proofs, and numerical results.

## Requirements & Environment
The analysis was performed in R and requires the following libraries:
```R
library(ismev)   # Extreme value modeling
library(evd)     # Extreme value distributions functions
library(tseries) # For ADF Test
