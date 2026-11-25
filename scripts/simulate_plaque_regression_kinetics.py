#!/usr/bin/env python3
"""
Multi-Pathway Plaque Regression & Dissolution Kinetic Simulation:
Comparing 5 Interventions over 24 Months (730 Days):
1. Untreated Control (Natural History)
2. Dietary Vinegar (Acetate, 15 mL/day)
3. EDTA Chelation Therapy (TACT/TACT2 Protocol)
4. Cholesterol-Crystal Solubilizer: 2HP-beta-Cyclodextrin (CD)
5. Ultra-Intensive Lipid Lowering (PCSK9i + High-Intensity Statin, GLAGOV/HUYGENS)
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

RESULTS_DIR = "/run/media/julian/ML6/arteriosklerosis/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

def plaque_dynamics_ode(y, t, params):
    """
    y[0]: L - Necrotic Lipid Core Volume (mm^3, baseline 50)
    y[1]: C - Calcification Volume (mm^3, baseline 20)
    y[2]: F - Fibrous Cap / Matrix Volume (mm^3, baseline 30)
    y[3]: CapThick - Fibrous Cap Thickness (um, baseline 55 um)
    y[4]: Inflam - Plaque Macrophage / NLRP3 / MMP-9 Activity (rel, baseline 1.0)
    """
    L, C, F, CapThick, Inflam = y
    
    # Kinetic parameters
    # Influx of new ApoB lipids into plaque depends on systemic LDL-C
    J_influx = params["J_influx"]        # mm^3/day
    k_efflux_base = params["k_efflux"]   # 1/day
    k_cd_solubilize = params["k_cd"]     # 1/day (direct cyclodextrin crystal dissolution)
    
    # Calcium calcification dynamics
    k_calc_growth = params["k_calc_growth"] # 1/day
    k_edta_dissolv = params["k_edta"]       # 1/day
    
    # Fibrous cap synthesis vs MMP degradation
    k_collagen_synth = params["k_synth"]    # um/day
    k_mmp_degrade = params["k_degrade"]     # um/day
    
    # 1. Lipid core rate of change
    # Net change = Influx - (Basal Efflux + CD-mediated Solubilization)
    dL_dt = J_influx - (k_efflux_base + k_cd_solubilize) * L
    
    # 2. Calcification rate of change
    # Proportional to local inflammation/osteogenic transdifferentiation minus EDTA chelation
    dC_dt = k_calc_growth * (0.5 + 0.5 * Inflam) - k_edta_dissolv * C
    
    # 3. Fibrous cap matrix volume
    dF_dt = 0.05 * (k_collagen_synth - k_mmp_degrade * Inflam)
    
    # 4. Fibrous cap thickness (um)
    # Target thickness balances VSMC collagen production and macrophage MMP digestion
    dCapThick_dt = k_collagen_synth - k_mmp_degrade * Inflam
    
    # 5. Inflammation dynamics
    # Stimulated by cholesterol crystals (L) and suppressed by interventions (AMPK, FFAR2, LXR)
    target_inflam = params["target_inflam_factor"] * (L / 50.0)
    dInflam_dt = 0.05 * (target_inflam - Inflam)
    
    return [dL_dt, dC_dt, dF_dt, dCapThick_dt, dInflam_dt]

# Parameter sets for each intervention
interventions = {
    "Control (Untreated)": {
        "J_influx": 0.045,           # High LDL influx (~130 mg/dL)
        "k_efflux": 0.0006,          # Weak natural baseline efflux
        "k_cd": 0.0,
        "k_calc_growth": 0.008,      # Progressive spotty microcalcification
        "k_edta": 0.0,
        "k_synth": 0.04,
        "k_degrade": 0.05,           # Active MMPs thin the cap
        "target_inflam_factor": 1.0,
        "color": "#d62728",          # Red
        "style": "-"
    },
    "Vinegar (15 mL/day)": {
        "J_influx": 0.042,           # Very slight reduction via transient postprandial glucose damping
        "k_efflux": 0.00065,         # Minor AMPK-mediated TFEB lipophagy
        "k_cd": 0.0,
        "k_calc_growth": 0.0078,     # Buffered in blood -> no direct calcium dissolution
        "k_edta": 0.0,
        "k_synth": 0.041,
        "k_degrade": 0.047,          # Mild endothelial NO anti-inflammatory effect
        "target_inflam_factor": 0.95,
        "color": "#ff7f0e",          # Orange
        "style": "--"
    },
    "EDTA Chelation Therapy": {
        "J_influx": 0.045,           # Does not lower LDL-C
        "k_efflux": 0.0006,
        "k_cd": 0.0,
        "k_calc_growth": 0.008,
        "k_edta": 0.0004,            # Minimal surface calcium extraction (TACT2 reality)
        "k_synth": 0.04,
        "k_degrade": 0.049,
        "target_inflam_factor": 0.98,
        "color": "#9467bd",          # Purple
        "style": "-."
    },
    "2HP-beta-Cyclodextrin": {
        "J_influx": 0.040,
        "k_efflux": 0.0012,          # LXR activation massively stimulates ABCA1/G1
        "k_cd": 0.0018,              # DIRECT cholesterol crystal solubilization
        "k_calc_growth": 0.004,
        "k_edta": 0.0,
        "k_synth": 0.05,
        "k_degrade": 0.025,          # Strong suppression of crystal-induced NLRP3
        "target_inflam_factor": 0.45,
        "color": "#1f77b4",          # Blue
        "style": "-"
    },
    "PCSK9i + Statin (Gold Standard)": {
        "J_influx": 0.004,           # LDL-C dropped to <30 mg/dL: influx practically halted!
        "k_efflux": 0.0015,          # Efflux strongly exceeds influx
        "k_cd": 0.0,
        "k_calc_growth": 0.009,      # Statins promote plaque macrocalcification (stabilization!)
        "k_edta": 0.0,
        "k_synth": 0.08,             # Thickening fibrous cap (+60 um, HUYGENS trial)
        "k_degrade": 0.015,          # Quenched vascular inflammation
        "target_inflam_factor": 0.35,
        "color": "#2ca02c",          # Green
        "style": "-"
    }
}

# Initial plaque conditions (Established vulnerable plaque)
# Total volume = 100 mm^3 (L=50, C=20, F=30)
# Cap thickness = 55 um (TCFA threshold <65 um)
# Inflammation = 1.0
y0 = [50.0, 20.0, 30.0, 55.0, 1.0]

# Time grid: 0 to 730 days (24 months)
t_days = np.linspace(0, 730, 730)
t_months = t_days / 30.417

solutions = {}
for name, p in interventions.items():
    sol = odeint(plaque_dynamics_ode, y0, t_days, args=(p,))
    total_vol = sol[:, 0] + sol[:, 1] + sol[:, 2]
    percent_tav_change = ((total_vol - total_vol[0]) / total_vol[0]) * 100.0
    solutions[name] = {
        "sol": sol,
        "total_vol": total_vol,
        "percent_tav_change": percent_tav_change
    }

# ==============================================================================
# PUBLICATION-GRADE DIAGNOSTIC FIGURE
# ==============================================================================
fig, axs = plt.subplots(2, 2, figsize=(16, 12), dpi=300)

# Panel A: Percent Change in Total Atheroma Volume (PAV / TAV)
ax1 = axs[0, 0]
for name, data in solutions.items():
    p = interventions[name]
    ax1.plot(t_months, data["percent_tav_change"], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax1.axhline(0, color='black', linestyle=':', lw=1.2)
ax1.set_title("A. Total Atheroma Volume (TAV) % Change", fontsize=12, fontweight='bold')
ax1.set_xlabel("Time (Months)", fontsize=11)
ax1.set_ylabel("Change in Plaque Volume (%)", fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9, loc='upper left')

# Panel B: Necrotic Lipid Core Volume (mm^3) - Actual Lipid Dissolution
ax2 = axs[0, 1]
for name, data in solutions.items():
    p = interventions[name]
    ax2.plot(t_months, data["sol"][:, 0], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax2.set_title("B. Necrotic Lipid Core Volume (Lipid Pool Regression)", fontsize=12, fontweight='bold')
ax2.set_xlabel("Time (Months)", fontsize=11)
ax2.set_ylabel("Lipid Core Volume (mm³)", fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=9, loc='upper right')

# Panel C: Fibrous Cap Thickness (um) - Plaque Stabilization Index
ax3 = axs[1, 0]
for name, data in solutions.items():
    p = interventions[name]
    ax3.plot(t_months, data["sol"][:, 3], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax3.axhline(65.0, color='red', linestyle='--', lw=1.8, label="Vulnerable Cap Threshold (<65 µm)")
ax3.set_title("C. Fibrous Cap Thickness (Rupture Resistance)", fontsize=12, fontweight='bold')
ax3.set_xlabel("Time (Months)", fontsize=11)
ax3.set_ylabel("Cap Thickness (µm)", fontsize=11)
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=9, loc='lower right')

# Panel D: Calcification Volume (mm^3) - Calcium Response
ax4 = axs[1, 1]
for name, data in solutions.items():
    p = interventions[name]
    ax4.plot(t_months, data["sol"][:, 1], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax4.set_title("D. Plaque Calcification Volume (Macro vs. Micro-calcification)", fontsize=12, fontweight='bold')
ax4.set_xlabel("Time (Months)", fontsize=11)
ax4.set_ylabel("Calcification Volume (mm³)", fontsize=11)
ax4.grid(True, alpha=0.3)
ax4.legend(fontsize=9, loc='upper left')

plt.tight_layout()
fig_path = os.path.join(RESULTS_DIR, "plaque_dissolution_comparison_figures.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"Generated comparison plot: {fig_path}")

# Generate numerical summary table
summary_results = {}
for name, data in solutions.items():
    sol = data["sol"]
    summary_results[name] = {
        "tav_change_24m_percent": round(float(data["percent_tav_change"][-1]), 2),
        "lipid_core_final_mm3": round(float(sol[-1, 0]), 2),
        "lipid_core_reduction_percent": round(float((50.0 - sol[-1, 0]) / 50.0 * 100.0), 2),
        "fibrous_cap_final_um": round(float(sol[-1, 3]), 1),
        "calcification_final_mm3": round(float(sol[-1, 1]), 2),
        "inflammation_final_rel": round(float(sol[-1, 4]), 2)
    }

with open(os.path.join(RESULTS_DIR, "plaque_dissolution_summary.json"), "w") as f:
    json.dump(summary_results, f, indent=2)
print("Saved plaque_dissolution_summary.json successfully.")
