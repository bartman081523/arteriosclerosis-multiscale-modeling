# 06. Glagov Remodeling, Luminal Coupling & Dual Reporting Framework

## 1. Mathematical Formalization of Luminal Feedback

### Glagov Outward Remodeling
In coronary arteries, early atheroma accumulation expands outward, preserving internal lumen area:

$$A_{\text{lumen}}(t) = A_{\text{EEM}}(t) - A_{\text{plaque}}(t)$$

- For Plaque Burden $PB \le 40\%$: External elastic lamina expands ($A_{\text{lumen}} \approx A_{\text{lumen},0}$).
- For Plaque Burden $PB > 40\%$: Outward reserve is exhausted; plaque encroaches directly on luminal cross-section.

### Poiseuille Shear Stress & Contact Inhibition
Luminal narrowing accelerates flow velocity, increasing wall shear stress:

$$\tau_w(t) = \tau_{w,0} \cdot \left(\frac{A_{\text{lumen},0}}{A_{\text{lumen}}(t)}\right)^{1.5}$$

Elevated WSS stimulates endothelial KLF2/NO, inhibiting VSMC synthetic overgrowth. Combined with the avascular oxygen diffusion limit ($T_{\text{cap,max}} \approx 120\ \mu\text{m}$), cap thickening is biochemically bounded:

$$v_{\text{synth}}(t) = v_{\text{VSMC},0} \cdot \frac{1}{1 + \left(\frac{\tau_w(t)}{\tau_{w,\text{threshold}}}\right)^2} \cdot \left(1 - \frac{T_{\text{cap}}(t)}{T_{\text{cap,max}}}\right)$$

---

## 2. Dual Reporting: Absolute vs. Relative Metrics (N=1000, 24 Months)

| Parameter | Control (Null) Median [95% CI] | Supermarket Protocol Median [95% CI] | Absolute Difference |
| :--- | :--- | :--- | :--- |
| **Plaque Volume (Absolute)** | **+16.98 mm³** [+4.64, +33.79] | **+0.63 mm³** [-7.92, +10.51] | **-16.35 mm³ avoided growth** |
| **Plaque Volume (Relative)** | **+17.32%** [+4.71%, +31.45%] | **+0.65%** [-7.90%, +10.23%] | **-16.67 percentage points** |
| **Fibrous Cap Thickness** | **54.7 µm** (stagnates in TCFA zone) | **68.1 µm** (stabilized above 65 µm) | **+13.4 µm thicker protective cap** |
| **Luminal Stenosis** | **0.0%** (Glagov compensated) | **0.0%** (Glagov compensated) | **No adverse lumen reduction** |

---

## 3. Re-Evaluated Sobol Total-Order Indices ($S_T$)
Under closed-loop luminal coupling:
- `k_deconj_macrophage`: $S_T = 0.799$
- `adherence_mean`: $S_T = 0.204$
- `vsmc_collagen_rate`: dropped from 1.000 to **0.954**
- `adherence_mean` on cap: **0.023** (active coupling)
- `mmp_degradation_rate` on cap: **0.017** (active coupling)
- `enos_hill_ec50` on cap: **0.007** (active coupling)

This confirms that the cap thickness is now fully embedded in the physiological network, responding dynamically to inflammation, adherence, and hemodynamics.
