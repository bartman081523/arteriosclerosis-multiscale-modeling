#!/usr/bin/env python3
"""
Multi-Scale Biophysical & PK/PD Modeling of Vinegar (Acetate) in Arteriosclerosis:
1. Blood Acid-Base & Plaque Dissolution Equilibrium (Henderson-Hasselbalch & Hydroxyapatite Ksp)
2. Pharmacokinetics of Oral Vinegar Ingestion (Compartment ODEs)
3. Endothelial & Macrophage FFAR2 Receptor Occupancy & eNOS/NO Dynamics
4. ACSS2 Intracellular Flux, AMP/ATP Shift, and AMPK Activation Kinetics
5. Downstream Atheroprotection: HMGCR Inactivation, Autophagy/TFEB, and Plaque Inflammation
6. Comparison: Oral Vinegar Bolus vs. Colonic Microbiome Fermentation
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

RESULTS_DIR = "/run/media/julian/ML6/arteriosklerosis/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ==============================================================================
# MODULE 1: ACID-BASE EQUILIBRIUM & CALCIFIED PLAQUE SOLUBILITY
# ==============================================================================
def simulate_acid_base_and_plaque(dose_ml=15.0):
    """
    Simulates systemic buffering and calcium phosphate solubility upon vinegar ingestion.
    """
    # 5% acetic acid solution (w/v): 50 g/L -> Molar mass = 60.05 g/mol -> 0.833 M
    acetic_acid_conc_M = 0.833
    dose_liters = dose_ml / 1000.0
    total_moles_acetate = dose_liters * acetic_acid_conc_M  # ~0.0125 mol (12.5 mmol for 15 ml)
    
    # Extracellular fluid (ECF) volume in adult human
    V_ecf = 14.0  # Liters
    
    # Basal physiological parameters
    pH_basal = 7.40
    pKa_acetate = 4.76
    pKa_bicarbonate = 6.10
    HCO3_basal = 0.024  # 24 mM (mol/L)
    dCO2_basal = 0.0012 # 1.2 mM (mol/L, pCO2 = 40 mmHg * 0.03)
    
    # Fraction of ionized acetate at pH 7.40
    # [Ac-] / [HAc] = 10^(7.40 - 4.76) = 10^2.64 = 436.5
    frac_ionized = 10**(pH_basal - pKa_acetate) / (1.0 + 10**(pH_basal - pKa_acetate))
    
    # Protons released into ECF if fully absorbed
    H_released_moles = total_moles_acetate * frac_ionized
    delta_HCO3 = H_released_moles / V_ecf
    
    # Immediate uncompensated pH
    HCO3_uncomp = max(0.005, HCO3_basal - delta_HCO3)
    pH_uncomp = pKa_bicarbonate + np.log10(HCO3_uncomp / dCO2_basal)
    delta_pH_uncomp = pH_uncomp - pH_basal
    
    # Compensated state (respiratory compensation hyperventilates CO2 out within minutes)
    # Winter's formula: expected pCO2 = 1.5 * [HCO3-] + 8
    pCO2_exp = 1.5 * (HCO3_uncomp * 1000) + 8
    dCO2_comp = pCO2_exp * 0.03 / 1000
    pH_comp = pKa_bicarbonate + np.log10(HCO3_uncomp / dCO2_comp)
    delta_pH_comp = pH_comp - pH_basal
    
    # Hydroxyapatite solubility product: Ca10(PO4)6(OH)2 <-> 10 Ca2+ + 6 PO4(3-) + 2 OH-
    # Ksp ~ 10^-116.8
    # Plaque dissolution occurs only when ion activity product IAP < Ksp (requires pH < 5.5)
    pH_range = np.linspace(5.0, 7.8, 100)
    # Relative solubility ratio S(pH) / S(7.4): log10(solubility) roughly increases ~1.5 units per pH drop below 6.5
    relative_solubility = 10**(1.5 * (7.4 - pH_range))
    
    return {
        "dose_ml": dose_ml,
        "total_mmol": total_moles_acetate * 1000,
        "frac_ionized_percent": frac_ionized * 100,
        "pH_basal": pH_basal,
        "pH_uncomp": pH_uncomp,
        "delta_pH_uncomp": delta_pH_uncomp,
        "pH_comp": pH_comp,
        "delta_pH_comp": delta_pH_comp,
        "pH_range": pH_range,
        "relative_solubility": relative_solubility
    }

# ==============================================================================
# MODULE 2: PHARMACOKINETICS & MULTI-COMPARTMENT DYNAMIC PATHWAY MODEL
# ==============================================================================
def pk_pd_ode_system(y, t, params):
    """
    Coupled ODE system:
    y[0]: Stomach acetate (mmol)
    y[1]: Portal vein / Gut acetate (mmol)
    y[2]: Systemic arterial plasma acetate concentration (uM)
    y[3]: Endothelial FFAR2 receptor activation (0 to 1)
    y[4]: Endothelial eNOS active fraction (phosphorylated Ser1177, normalized)
    y[5]: Vascular Nitric Oxide (NO) release rate (relative)
    y[6]: Macrophage FFAR2 activation (0 to 1)
    y[7]: Macrophage NLRP3 Inflammasome activity (normalized, basal = 1.0)
    y[8]: Cellular ACSS2 flux (umol/L/min)
    y[9]: Cellular AMP/ATP ratio (normalized, basal = 1.0)
    y[10]: AMPK activation state (phosphorylated Thr172, normalized, basal = 1.0)
    y[11]: HMG-CoA Reductase activity (normalized, basal = 1.0)
    y[12]: Macrophage TFEB Nuclear Autophagy score (normalized, basal = 1.0)
    """
    (A_stom, A_gut, C_plasma, FFAR2_endo, eNOS_act, NO_rel, 
     FFAR2_mac, NLRP3_act, ACSS2_flux, AMP_ATP, AMPK_act, HMGCR_act, TFEB_autoph) = y
    
    # Parameters
    k_empty = params["k_empty"]       # Gastric emptying rate (1/min)
    k_abs = params["k_abs"]           # Intestinal absorption rate (1/min)
    f_first_pass = params["f_first_pass"] # Hepatic first-pass extraction (0.65)
    V_plasma = params["V_plasma"]     # Volume of distribution (L, ~12 L ECF accessible)
    k_elim = params["k_elim"]         # Systemic metabolic elimination (1/min)
    C_basal = params["C_basal"]       # Basal endogenous plasma acetate (uM, ~75 uM)
    
    # Hill constants
    Kd_FFAR2 = params["Kd_FFAR2"]     # ~350 uM
    Km_ACSS2 = params["Km_ACSS2"]     # ~70 uM
    Vmax_ACSS2 = params["Vmax_ACSS2"] # ~15 umol/L/min
    
    # 1. Pharmacokinetics
    dA_stom_dt = - k_empty * A_stom
    dA_gut_dt = k_empty * A_stom - k_abs * A_gut
    
    # Systemic input: gut absorption minus hepatic first pass
    input_systemic_mmol_min = (1.0 - f_first_pass) * k_abs * A_gut
    # Convert mmol/min into uM/min in V_plasma: (mmol / L) * 1000 = uM
    dC_plasma_dt = (input_systemic_mmol_min / V_plasma) * 1000.0 - k_elim * (C_plasma - C_basal)
    
    # 2. FFAR2 Occupancy (quasi-steady state or rapid equilibrium)
    target_FFAR2 = C_plasma / (C_plasma + Kd_FFAR2)
    dFFAR2_endo_dt = 1.5 * (target_FFAR2 - FFAR2_endo)
    dFFAR2_mac_dt = 1.5 * (target_FFAR2 - FFAR2_mac)
    
    # 3. Endothelial eNOS & NO Signaling
    # Activation by FFAR2 (via Ca2+/CaMKK and AMPK)
    target_eNOS = 1.0 + 0.8 * FFAR2_endo + 0.4 * (AMPK_act - 1.0)
    deNOS_act_dt = 0.2 * (target_eNOS - eNOS_act)
    dNO_rel_dt = 0.5 * (eNOS_act - NO_rel)
    
    # 4. Macrophage NLRP3 Inflammasome
    # Suppressed by FFAR2 signaling and autophagy
    target_NLRP3 = 1.0 / (1.0 + 1.2 * FFAR2_mac + 0.5 * (TFEB_autoph - 1.0))
    dNLRP3_act_dt = 0.1 * (target_NLRP3 - NLRP3_act)
    
    # 5. Cellular ACSS2 & AMP/ATP Ratio
    target_ACSS2_flux = Vmax_ACSS2 * C_plasma / (C_plasma + Km_ACSS2)
    dACSS2_flux_dt = 0.5 * (target_ACSS2_flux - ACSS2_flux)
    
    # ACSS2 uses ATP -> generates AMP + PPi
    # Shift in AMP/ATP ratio is proportional to ACSS2 flux above basal
    basal_flux = Vmax_ACSS2 * C_basal / (C_basal + Km_ACSS2)
    target_AMP_ATP = 1.0 + 0.05 * max(0.0, ACSS2_flux - basal_flux)
    dAMP_ATP_dt = 0.2 * (target_AMP_ATP - AMP_ATP)
    
    # 6. AMPK Activation
    # AMPK is allosterically activated and phosphorylated by elevated AMP/ATP
    target_AMPK = 1.0 + 1.2 * (AMP_ATP - 1.0)
    dAMPK_act_dt = 0.15 * (target_AMPK - AMPK_act)
    
    # 7. HMG-CoA Reductase Inactivation (Phosphorylation by AMPK at Ser872)
    # Active HMGCR fraction decreases as AMPK increases
    target_HMGCR = 1.0 / (1.0 + 0.9 * (AMPK_act - 1.0))
    dHMGCR_act_dt = 0.15 * (target_HMGCR - HMGCR_act)
    
    # 8. TFEB Nuclear Translocation & Macrophage Autophagy
    # Phosphorylation of ACSS2 by AMPK promotes nuclear translocation and TFEB-mediated transcription
    target_TFEB = 1.0 + 0.7 * (AMPK_act - 1.0) * (ACSS2_flux / basal_flux)
    dTFEB_autoph_dt = 0.08 * (target_TFEB - TFEB_autoph)
    
    return [dA_stom_dt, dA_gut_dt, dC_plasma_dt, dFFAR2_endo_dt, deNOS_act_dt, dNO_rel_dt,
            dFFAR2_mac_dt, dNLRP3_act_dt, dACSS2_flux_dt, dAMP_ATP_dt, dAMPK_act_dt, dHMGCR_act_dt, dTFEB_autoph_dt]

def run_simulation():
    # Parameters
    params = {
        "k_empty": 0.04,          # Gastric emptying: ~25 min half-time
        "k_abs": 0.08,            # Rapid intestinal absorption
        "f_first_pass": 0.65,     # 65% hepatic extraction
        "V_plasma": 12.0,         # Liters distribution volume
        "k_elim": 0.06,           # Elimination rate: half-life ~12 min
        "C_basal": 75.0,          # 75 uM basal plasma acetate
        "Kd_FFAR2": 350.0,        # 350 uM EC50/Kd
        "Km_ACSS2": 70.0,         # 70 uM Km
        "Vmax_ACSS2": 15.0        # 15 umol/L/min
    }
    
    # Time grid: 0 to 240 minutes (4 hours)
    t = np.linspace(0, 240, 500)
    
    # Condition 1: 15 mL 5% vinegar (12.5 mmol acetate)
    y0_15ml = [12.5, 0.0, 75.0, 75.0/(75.0+350.0), 1.0, 1.0, 75.0/(75.0+350.0), 1.0, 
               15.0*75.0/(75.0+70.0), 1.0, 1.0, 1.0, 1.0]
    sol_15ml = odeint(pk_pd_ode_system, y0_15ml, t, args=(params,))
    
    # Condition 2: 30 mL 5% vinegar (25.0 mmol acetate)
    y0_30ml = [25.0, 0.0, 75.0, 75.0/(75.0+350.0), 1.0, 1.0, 75.0/(75.0+350.0), 1.0, 
               15.0*75.0/(75.0+70.0), 1.0, 1.0, 1.0, 1.0]
    sol_30ml = odeint(pk_pd_ode_system, y0_30ml, t, args=(params,))
    
    # Acid-Base calculations
    ab_15 = simulate_acid_base_and_plaque(15.0)
    ab_30 = simulate_acid_base_and_plaque(30.0)
    
    # ==============================================================================
    # PLOTTING PUBLICATION-GRADE DIAGNOSTICS
    # ==============================================================================
    fig = plt.figure(figsize=(18, 12), dpi=300)
    
    # Panel 1: Pharmacokinetics - Plasma Acetate Concentration (uM)
    ax1 = fig.add_subplot(2, 3, 1)
    ax1.plot(t, sol_15ml[:, 2], 'b-', lw=2.5, label="15 mL Vinegar (12.5 mmol)")
    ax1.plot(t, sol_30ml[:, 2], 'b--', lw=2.5, label="30 mL Vinegar (25.0 mmol)")
    ax1.axhline(params["C_basal"], color='gray', linestyle=':', label="Basal Microbiome Level (75 µM)")
    ax1.axhline(params["Kd_FFAR2"], color='red', linestyle='-.', alpha=0.7, label="FFAR2 EC50 (350 µM)")
    ax1.axhline(params["Km_ACSS2"], color='green', linestyle='-.', alpha=0.7, label="ACSS2 Km (70 µM)")
    ax1.set_title("A. Systemic Plasma Acetate Pharmacokinetics", fontsize=12, fontweight='bold')
    ax1.set_xlabel("Time (minutes)", fontsize=11)
    ax1.set_ylabel("Plasma Acetate [µM]", fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9, loc='upper right')
    
    # Panel 2: Receptor Occupancy & Endothelial eNOS / NO Release
    ax2 = fig.add_subplot(2, 3, 2)
    ax2.plot(t, sol_15ml[:, 3] * 100, 'purple', lw=2.0, label="FFAR2 Occupancy (%) [15 mL]")
    ax2.plot(t, sol_15ml[:, 5], 'g-', lw=2.5, label="Endothelial NO Release (rel.) [15 mL]")
    ax2.plot(t, sol_30ml[:, 5], 'g--', lw=2.5, label="Endothelial NO Release (rel.) [30 mL]")
    ax2.set_title("B. Endothelial FFAR2 & eNOS/NO Activation", fontsize=12, fontweight='bold')
    ax2.set_xlabel("Time (minutes)", fontsize=11)
    ax2.set_ylabel("Occupancy (%) / Relative NO", fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc='upper right')
    
    # Panel 3: ACSS2 Flux, AMP/ATP Ratio, and AMPK Activation
    ax3 = fig.add_subplot(2, 3, 3)
    ax3.plot(t, sol_15ml[:, 10], 'orange', lw=2.5, label="AMPK Activity (Thr172-P) [15 mL]")
    ax3.plot(t, sol_30ml[:, 10], 'darkorange', lw=2.5, linestyle='--', label="AMPK Activity [30 mL]")
    ax3.plot(t, sol_15ml[:, 9], 'brown', lw=1.8, linestyle=':', label="AMP/ATP Ratio [15 mL]")
    ax3.set_title("C. ACSS2-Mediated AMPK Activation Cascade", fontsize=12, fontweight='bold')
    ax3.set_xlabel("Time (minutes)", fontsize=11)
    ax3.set_ylabel("Normalized Ratio / Activity", fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=9, loc='upper right')
    
    # Panel 4: Downstream Anti-Atherogenic Actions (HMGCR Inactivation & NLRP3 Suppression)
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.plot(t, sol_15ml[:, 11], 'r-', lw=2.5, label="HMG-CoA Reductase Activity (15 mL)")
    ax4.plot(t, sol_15ml[:, 7], 'm-', lw=2.5, label="NLRP3 Inflammasome Activity (15 mL)")
    ax4.plot(t, sol_15ml[:, 12], 'c-', lw=2.5, label="Macrophage TFEB Autophagy (15 mL)")
    ax4.axhline(1.0, color='gray', linestyle=':')
    ax4.set_title("D. Downstream Cellular Protection Targets", fontsize=12, fontweight='bold')
    ax4.set_xlabel("Time (minutes)", fontsize=11)
    ax4.set_ylabel("Activity relative to Basal (1.0)", fontsize=11)
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=9, loc='center right')
    
    # Panel 5: Biophysical Acid-Base & Plaque Solubility Falsification
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.plot(ab_15["pH_range"], ab_15["relative_solubility"], 'k-', lw=2.5, label="Hydroxyapatite Solubility")
    ax5.axvline(7.40, color='blue', linestyle='--', lw=2.0, label="Physiological pH (7.40)")
    ax5.axvline(ab_15["pH_uncomp"], color='magenta', linestyle=':', lw=2.0, 
                label=f"Vinegar 15ml Uncomp (pH {ab_15['pH_uncomp']:.2f})")
    ax5.axvspan(5.0, 5.8, color='red', alpha=0.2, label="Required Plaque Dissolution Zone (<5.8)")
    ax5.set_yscale("log")
    ax5.set_title("E. Falsification: Plaque Dissolution vs. pH", fontsize=12, fontweight='bold')
    ax5.set_xlabel("Blood / Tissue pH", fontsize=11)
    ax5.set_ylabel("Relative Hydroxyapatite Solubility (log)", fontsize=11)
    ax5.grid(True, alpha=0.3)
    ax5.legend(fontsize=8, loc='upper right')
    
    # Panel 6: Dietary Vinegar vs. Gut Microbiota Fermentation Daily Yield
    ax6 = fig.add_subplot(2, 3, 6)
    categories = ['15 mL Vinegar\n(1 tbsp)', '30 mL Vinegar\n(2 tbsp)', 'Colonic Microbiota\n(Low Fiber Diet)', 'Colonic Microbiota\n(High Fiber Diet)']
    daily_acetate_mmol = [12.5, 25.0, 120.0, 280.0]
    colors = ['skyblue', 'steelblue', 'goldenrod', 'darkgreen']
    bars = ax6.bar(categories, daily_acetate_mmol, color=colors, edgecolor='black', width=0.6)
    for bar in bars:
        yval = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2.0, yval + 5.0, f"{yval:.1f} mmol", ha='center', va='bottom', fontweight='bold', fontsize=9)
    ax6.set_title("F. Daily Acetate Source Comparison", fontsize=12, fontweight='bold')
    ax6.set_ylabel("Daily Acetate Production (mmol/day)", fontsize=11)
    ax6.set_ylim(0, 320)
    ax6.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plot_path = os.path.join(RESULTS_DIR, "arteriosclerosis_modeling_figures.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Generated multi-panel publication plot: {plot_path}")
    
    # Return numerical summary dictionary
    summary = {
        "c_max_15ml_uM": float(np.max(sol_15ml[:, 2])),
        "t_max_15ml_min": float(t[np.argmax(sol_15ml[:, 2])]),
        "c_max_30ml_uM": float(np.max(sol_30ml[:, 2])),
        "t_max_30ml_min": float(t[np.argmax(sol_30ml[:, 2])]),
        "max_ffar2_occupancy_15ml": float(np.max(sol_15ml[:, 3])),
        "max_ffar2_occupancy_30ml": float(np.max(sol_30ml[:, 3])),
        "max_no_increase_15ml_percent": float((np.max(sol_15ml[:, 5]) - 1.0) * 100),
        "max_hmgcr_inhibition_15ml_percent": float((1.0 - np.min(sol_15ml[:, 11])) * 100),
        "max_nlrp3_suppression_15ml_percent": float((1.0 - np.min(sol_15ml[:, 7])) * 100),
        "max_tfeb_autophagy_increase_15ml_percent": float((np.max(sol_15ml[:, 12]) - 1.0) * 100),
        "acid_base": {
            "uncompensated_delta_pH_15ml": ab_15["delta_pH_uncomp"],
            "compensated_delta_pH_15ml": ab_15["delta_pH_comp"],
            "uncompensated_delta_pH_30ml": ab_30["delta_pH_uncomp"],
            "compensated_delta_pH_30ml": ab_30["delta_pH_comp"],
        }
    }
    
    import json
    with open(os.path.join(RESULTS_DIR, "simulation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("Saved simulation_summary.json successfully.")

if __name__ == "__main__":
    run_simulation()
