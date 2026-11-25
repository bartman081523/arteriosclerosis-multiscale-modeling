# 05. The 4-Stage Verification Gate & Global Sobol Sensitivity Analysis

## 1. Verification Protocol (The 4 Gates)
To expose hidden dependencies, numerical artifacts, and over-precision, the simulation engine was tested through a formal 4-stage audit:

| Gate | Type | Method | Acceptance Criteria | Result |
| :--- | :--- | :--- | :--- | :--- |
| **1. Numerics** | Convergence | $N \in \{250, 1000, 5000\}$, $\Delta t \in \{0.5, 1.0\}\text{ d}$, 5 seeds | $\Delta_{\text{step}} < 0.30\%$, $\sigma_{\text{seed}} < 0.35\%$ | **PASSED** ($\Delta_{\text{dt}} = 0.000\%$, $\sigma_{\text{seed}} = 0.227\%$) |
| **2. Falsification** | Null Model | Placebo mode ($active=0$) under full noise | Matches untreated control; no NaNs | **PASSED** (Plaque: +17.12%, bounds intact) |
| **3. Ablation** | LOMO | Leave-One-Mechanism-Out | Quantify individual module contributions | **Identified Dominance** of $k_{\text{deconj}}$ |
| **4. Sensitivity** | Sobol | Saltelli sample ($N=256$, 3,584 evaluations) | Decompose variance into $S_1$, $S_T$ | **Exposed $S_1 = 0.775$** for deconjugation |

---

## 2. Submodel Ablation (LOMO) Results
- **Full Model:** Plaque median = **+0.65%**, Cap = **90.3 µm**.
- **Without Macrophage Deconjugation:** Plaque median jumps to **+7.20% ($\Delta = +6.55\%$)**.
- **Without SREBP-2 Feedback:** Plaque median = **+0.07% ($\Delta = -0.58\%$)**.
- **Without Hill Saturation:** Plaque median = **+0.65% ($\Delta = 0.00\%$)**.
- **Without Infection Spikes:** Plaque median = **+0.65% ($\Delta = 0.00\%$)**.

---

## 3. Global Sobol Sensitivity Indices (SALib)
Across 3,584 model evaluations:
- **Plaque Volume Median Variance:**
  - `k_deconj_macrophage`: $S_1 = 0.775, S_T = 0.771$ (Explains **77.5% of total variance**!).
  - `adherence_mean`: $S_1 = 0.191, S_T = 0.212$ (Explains **19.1% of variance**).
  - Together, local deconjugation and adherence account for **96.6%** of all plaque volume outcomes.
- **Interquartile Range (IQR / Heterogeneity):**
  - `k_deconj_macrophage`: $S_T = 0.970$ (**97.0% of between-patient spread** is driven by deconjugation efficiency!).
- **Cap Thickness Variance (Uncoupled Model):**
  - `vsmc_collagen_rate`: $S_1 = 0.997, S_T = 1.000$ (Revealed an uncoupled 1D artifact requiring luminal feedback).
