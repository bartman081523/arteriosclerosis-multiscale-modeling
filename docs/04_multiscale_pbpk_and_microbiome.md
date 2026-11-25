# 04. Multiscale PBPK & Microbiome-Stratified Population Architecture

## 1. PBPK Modeling & Phase-II Conjugation
In vivo, dietary polyphenols undergo rapid first-pass enteral and hepatic glucuronidation (UGTs) and sulfation (SULTs):
- **Free Hydroxytyrosol:** Circulates at only **20–50 nM** ($t_{1/2} \approx 45\ \text{min}$).
- **Conjugated Metabolites:** Reach **1.5–3.2 µM**.
- **Plaque Microenvironment Deconjugation:** In active atheromas, activated macrophages release lysosomal **$\beta$-glucuronidase** and **arylsulfatase**, partially re-activating free catechols locally at the inflammatory front.

---

## 2. Gut Microbiome Stratification
Metabolic conversion of pomegranate ellagitannins depends on colonic colonization by *Gordonibacter* and *Ellagibacter*:
- **Metabotype A (36.6%):** Efficient Urolithin A producers ($C_{\text{plasma}} \approx 0.5–1.8\ \mu\text{M}$) $\rightarrow$ full PON1 and autophagy stimulation.
- **Metabotype B (45.9%):** Produce Isourolithin A and Urolithin B.
- **Metabotype 0 (17.5%):** Non-producers (undetectable serum Urolithin A).

---

## 3. Homeostatic Counter-Regulation & Saturation
- **eNOS Hill Capping:** Combined agonism by acetate, SAC, and hydroxytyrosol saturates at maximal enzymatic capacity ($V_{\text{max}} \approx 1.7\times$), preventing unphysical linear summation.
- **Hepatic SREBP-2 Compensation:** Bile acid sequestration by $\beta$-glucan depletes intrahepatic cholesterol, triggering SCAP/SREBP-2 nuclear translocation and elevating compensatory HMG-CoA reductase synthesis by **+18–25%**, dampening net LDL reduction from theoretical -22% to realistic **-8% to -11%**.

---

## 4. Virtual Patient Population (N=1000, 730 Days)
- **Stochastic Daily Adherence:** Modeled via a Markov/Bernoulli process with individual adherence $p_{\text{adh}} \sim \mathcal{U}(0.68, 0.88)$ (cohort mean: **77.9%**).
- **Acute Inflammatory Perturbations:** Stochastically generated viral/bacterial infections (Poisson rate $\lambda = 2.0/\text{year}$) producing 5-day hs-CRP spikes (+6 to +12 mg/L) and transient MMP surges.
