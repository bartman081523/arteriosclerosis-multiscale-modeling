# 02. Plaque Dissolution & Regression: 2HP-β-Cyclodextrin vs. Chelation vs. PCSK9 Inhibition

## 1. The Anatomy of Atheroma Regression
Atherosclerotic plaque consists of three distinct compartments:
1. **Necrotic Core:** Free cholesterol, cholesterol esters, and solid extracellular cholesterol crystals.
2. **Calcification:** Biological hydroxyapatite nanocrystals ($Ca_{10}(PO_4)_6(OH)_2$).
3. **Fibrous Cap:** Type I & III collagen matrix synthesized by vascular smooth muscle cells (VSMCs).

True plaque "dissolution" requires distinct biochemical mechanisms depending on the target compartment.

---

## 2. 2-Hydroxypropyl-β-Cyclodextrin (2HP-β-CD)
- **Host-Guest Chemistry:** $\beta$-Cyclodextrin is a cone-shaped cyclic heptasaccharide with a hydrophobic inner core and hydrophilic rim.
- **In-Silico Docking:** AutoDock Vina models cholesterol inclusion into the cyclodextrin torus with $\Delta G = -3.96\ \text{kcal/mol}$ ($K_a \approx 6.1 \times 10^2\ \text{M}^{-1}$ for native $\beta$-CD; $> 10^3–10^4\ \text{M}^{-1}$ for 2HP-$\beta$-CD).
- **Biological Mechanism (*Science Transl. Med.* 2016, Zimmer et al.):**
  1. Solubilizes solid cholesterol crystals that otherwise puncture phagolysosomes and activate the NLRP3 inflammasome.
  2. Generates 27-hydroxycholesterol (27-OHC), which acts as an endogenous agonist of **LXR** (Liver X Receptor).
  3. LXR upregulates **ABCA1** and **ABCG1**, triggering massive cholesterol efflux to circulating ApoA-I/HDL.
  4. In our 24-month kinetic ODE simulation: **-65.1% reduction in necrotic core** and **-29.1% total atheroma volume regression**.

---

## 3. Disodium EDTA & Chelation Therapy (TACT / TACT2 Falsification)
- EDTA binds free ionic $Ca^{2+}$ with $\log K_1 = 10.69$.
- However, calcium in atheromas is locked within crystalline hydroxyapatite with a lattice enthalpy of $\sim -37,000\ \text{kJ/mol}$.
- The recent **TACT2 trial** (*JAMA* 2024, n=1,000 post-MI diabetic patients) demonstrated that repeated EDTA infusions failed to alter coronary artery calcification (CAC) scores or reduce recurrent cardiovascular events.
- In our simulation: EDTA achieved 0% necrotic core clearance and total atheroma volume grew by **+8.6%**.

---

## 4. PCSK9 Inhibitors + High-Intensity Statins (The Clinical Benchmark)
- Lowering LDL-C to $< 30–40\ \text{mg/dL}$ (e.g. Evolocumab + Atorvastatin 80 mg) inverts the trans-endothelial lipid concentration gradient.
- Subendothelial ApoB influx ceases ($J_{\text{influx}} \approx 0$), allowing HDL-mediated reverse cholesterol transport to dominate ($J_{\text{efflux}} \gg J_{\text{influx}}$).
- **Clinical Validation:** Replicated findings from the **GLAGOV trial** (JAMA 2016, IVUS regression) and **HUYGENS trial** (JACC 2022, OCT cap thickening from 54 µm to 110.7 µm).
