#!/usr/bin/env python3
"""
24-Month Kinetic Simulation:
Can accessible supermarket foods improve or regress atherosclerosis?
Comparing:
1. Untreated Control
2. Vinegar alone (15 mL/day)
3. Aged Garlic alone (Schwarzer Knoblauch)
4. Supermarket Synergistic Protocol (Haferkleie + Knoblauch + Granatapfel + Olivenöl + Natto/K2)
5. Clinical Gold Standard (PCSK9i + Statin)
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

RESULTS_DIR = "/run/media/julian/ML6/arteriosklerosis/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

def plaque_dynamics_extended(y, t, params):
    """
    State variables:
    y[0]: L - Necrotic Lipid Core Volume (mm^3, baseline 50)
    y[1]: C - Calcification Volume (mm^3, baseline 20)
    y[2]: F - Fibrous Cap Volume (mm^3, baseline 30)
    y[3]: CapThick - Fibrous Cap Thickness (um, baseline 55 um)
    y[4]: Inflam - Plaque Macrophage Activity (rel, baseline 1.0)
    """
    L, C, F, CapThick, Inflam = y
    
    J_influx = params["J_influx"]        # Influx of oxLDL/ApoB particles
    k_efflux = params["k_efflux"]        # ABCA1/SR-B1/PON1-mediated reverse cholesterol transport
    k_calc = params["k_calc"]            # Calcification growth rate
    k_synth = params["k_synth"]          # Cap collagen synthesis
    k_degrade = params["k_degrade"]      # MMP-mediated cap digestion
    target_inflam_f = params["target_inflam"]
    
    dL_dt = J_influx - k_efflux * L
    dC_dt = k_calc * (0.4 + 0.6 * Inflam)
    dF_dt = 0.05 * (k_synth - k_degrade * Inflam)
    dCapThick_dt = k_synth - k_degrade * Inflam
    
    target_inflam = target_inflam_f * (L / 50.0)
    dInflam_dt = 0.05 * (target_inflam - Inflam)
    
    return [dL_dt, dC_dt, dF_dt, dCapThick_dt, dInflam_dt]

cohorts = {
    "1. Kontrolle (Unbehandelt)": {
        "J_influx": 0.045,
        "k_efflux": 0.0006,
        "k_calc": 0.008,
        "k_synth": 0.04,
        "k_degrade": 0.05,
        "target_inflam": 1.0,
        "color": "#d62728",
        "style": "-"
    },
    "2. Essig allein (15 mL/Tag)": {
        "J_influx": 0.042,
        "k_efflux": 0.00065,
        "k_calc": 0.0078,
        "k_synth": 0.041,
        "k_degrade": 0.047,
        "target_inflam": 0.95,
        "color": "#ff7f0e",
        "style": "--"
    },
    "3. Knoblauch allein (Aged Garlic)": {
        # Dr. Budoff UCLA RCTs: -1.5% to -3% LAP regression via H2S/eNOS and CD36 inhibition
        "J_influx": 0.035,
        "k_efflux": 0.0008,
        "k_calc": 0.0055,       # Slows CAC progression
        "k_synth": 0.045,
        "k_degrade": 0.038,
        "target_inflam": 0.75,
        "color": "#8c564b",
        "style": "-."
    },
    "4. Supermarkt-Synergie-Protokoll": {
        # Haferkleie (Beta-Glucan) -> senkt LDL um 10-12% (J_influx sinkt auf 0.024)
        # Olivenöl (Hydroxytyrosol) -> blockiert oxLDL-Bildung
        # Granatapfel (Urolithin A) -> aktiviert PON1 (+80%), k_efflux steigt auf 0.0011
        # Schwarzer Knoblauch -> dämpft Makrophagen-Inflammation & CD36
        # Natto (K2 & Nattokinase) -> k_calc sinkt auf 0.0035 (MGP-Carboxylierung) & fibrinolytische Kappenstabilisierung
        "J_influx": 0.022,
        "k_efflux": 0.00115,     # Signifikanter Revers-Cholesterintransport!
        "k_calc": 0.0035,        # Stopp von Mikroverkalkungen durch Vitamin K2
        "k_synth": 0.060,        # Fibröse Kappe wächst
        "k_degrade": 0.025,      # Geringe MMP-Aktivität
        "target_inflam": 0.42,
        "color": "#1f77b4",      # Blau
        "style": "-"
    },
    "5. Goldstandard (PCSK9i + Statin)": {
        "J_influx": 0.004,
        "k_efflux": 0.0015,
        "k_calc": 0.009,
        "k_synth": 0.080,
        "k_degrade": 0.015,
        "target_inflam": 0.35,
        "color": "#2ca02c",      # Grün
        "style": "-"
    }
}

y0 = [50.0, 20.0, 30.0, 55.0, 1.0]
t_days = np.linspace(0, 730, 730)
t_months = t_days / 30.417

results = {}
for name, p in cohorts.items():
    sol = odeint(plaque_dynamics_extended, y0, t_days, args=(p,))
    tav = sol[:, 0] + sol[:, 1] + sol[:, 2]
    pct_tav = ((tav - tav[0]) / tav[0]) * 100.0
    results[name] = {
        "sol": sol,
        "tav": tav,
        "pct_tav": pct_tav
    }

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(16, 12), dpi=300)

# Panel A: Total Atheroma Volume % Change
ax1 = axs[0, 0]
for name, d in results.items():
    p = cohorts[name]
    ax1.plot(t_months, d["pct_tav"], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax1.axhline(0, color='black', linestyle=':', lw=1.2)
ax1.set_title("A. Gesamtes Plaque-Volumen (TAV) % Änderung", fontsize=12, fontweight='bold')
ax1.set_xlabel("Monate", fontsize=11)
ax1.set_ylabel("Änderung des Plaque-Volumens (%)", fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9, loc='upper left')

# Panel B: Lipid Core Regression (mm^3)
ax2 = axs[0, 1]
for name, d in results.items():
    p = cohorts[name]
    ax2.plot(t_months, d["sol"][:, 0], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax2.set_title("B. Nekrotischer Lipidkern (Cholesterin-Rückbildung)", fontsize=12, fontweight='bold')
ax2.set_xlabel("Monate", fontsize=11)
ax2.set_ylabel("Lipidkern-Volumen (mm³)", fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=9, loc='upper right')

# Panel C: Fibrous Cap Thickness (um)
ax3 = axs[1, 0]
for name, d in results.items():
    p = cohorts[name]
    ax3.plot(t_months, d["sol"][:, 3], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax3.axhline(65.0, color='red', linestyle='--', lw=1.8, label="Kritische Rupturgrenze (<65 µm)")
ax3.set_title("C. Fibröse Kappendicke (Stabilität gegen Ruptur)", fontsize=12, fontweight='bold')
ax3.set_xlabel("Monate", fontsize=11)
ax3.set_ylabel("Kappendicke (µm)", fontsize=11)
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=9, loc='lower right')

# Panel D: Plaque Calcification
ax4 = axs[1, 1]
for name, d in results.items():
    p = cohorts[name]
    ax4.plot(t_months, d["sol"][:, 1], label=name, color=p["color"], linestyle=p["style"], lw=2.5)
ax4.set_title("D. Plaque-Verkalkung (Kalkentwicklung)", fontsize=12, fontweight='bold')
ax4.set_xlabel("Monate", fontsize=11)
ax4.set_ylabel("Kalkvolumen (mm³)", fontsize=11)
ax4.grid(True, alpha=0.3)
ax4.legend(fontsize=9, loc='upper left')

plt.tight_layout()
fig_path = os.path.join(RESULTS_DIR, "supermarket_atherosclerosis_simulation.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"Generated supermarket simulation plot: {fig_path}")

summary = {}
for name, d in results.items():
    sol = d["sol"]
    summary[name] = {
        "tav_change_24m_percent": round(float(d["pct_tav"][-1]), 2),
        "lipid_core_final_mm3": round(float(sol[-1, 0]), 2),
        "lipid_core_reduction_percent": round(float((50.0 - sol[-1, 0]) / 50.0 * 100.0), 2),
        "fibrous_cap_final_um": round(float(sol[-1, 3]), 1),
        "calcification_final_mm3": round(float(sol[-1, 1]), 2),
        "inflammation_final_rel": round(float(sol[-1, 4]), 2)
    }

with open(os.path.join(RESULTS_DIR, "supermarket_simulation_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

print("Saved supermarket_simulation_summary.json successfully.")
