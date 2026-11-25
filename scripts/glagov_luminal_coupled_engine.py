#!/usr/bin/env python3
"""
Coupled Glagov-Remodeling, Luminal Shear Feedback, and Multi-Scale Plaque Engine:
- Glagov Law: Compensatory outward remodeling up to 40% plaque burden
- Luminal Narrowing & Poiseuille Shear Stress: tau_w(t) = tau_w0 * (A_lumen0 / A_lumen(t))^1.5
- Mechanosensitive feedback & Cap Saturation: Elevated WSS and O2 diffusion limit cap growth
- Dual reporting: Absolute shifts (mm^3, um, %-points) alongside relative rates
- Sobol sensitivity re-evaluation on coupled T_cap
"""
import os
import json
from dataclasses import dataclass
from typing import Dict, Tuple
import numpy as np
import matplotlib.pyplot as plt

RESULTS_DIR = "/run/media/julian/ML6/arteriosklerosis/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

@dataclass(frozen=True)
class CoupledConfig:
    n_patients: int = 1000
    days: int = 730
    dt: float = 1.0
    seed: int = 42
    placebo_null_mode: bool = False

@dataclass(frozen=True)
class CoupledParams:
    k_deconj_macrophage: float = 0.35
    srebp2_rebound_gain: float = 0.22
    enos_hill_ec50: float = 1.20
    vsmc_collagen_rate: float = 0.045
    mmp_degradation_rate: float = 0.038
    adherence_mean: float = 0.78
    # Glagov & Luminal Coupling Parameters
    glagov_threshold: float = 0.40      # 40% plaque burden outward remodeling limit
    cap_max_diffusion: float = 120.0    # um, max cap thickness before hypoxia limits collagen
    wss_inhibition_ec50: float = 2.0    # Pa, WSS threshold inhibiting VSMC overgrowth

def simulate_coupled_cohort(cfg: CoupledConfig, p: CoupledParams) -> Dict[str, np.ndarray]:
    rng = np.random.default_rng(cfg.seed)
    n = cfg.n_patients
    steps = int(cfg.days / cfg.dt)

    # Patient baseline
    metabotype = rng.choice([0, 1, 2], size=n, p=[0.366, 0.459, 0.175])
    urolithin_factor = np.select([metabotype == 0, metabotype == 1, metabotype == 2], [1.0, 0.45, 0.0])
    adherence = np.clip(rng.normal(p.adherence_mean, 0.06, size=n), 0.50, 1.0)
    base_ldl = np.clip(rng.normal(138.0, 24.0, size=n), 90.0, 210.0)
    base_crp = np.clip(rng.lognormal(0.8, 0.4, size=n), 0.5, 8.0)
    
    # Vessel geometry: baseline coronary artery (LAD/RCA)
    # Vessel total area A_EEM0 ~ 15 mm^2, baseline lumen ~ 10 mm^2, baseline plaque ~ 5 mm^2 (burden ~ 33%)
    a_eem_0 = np.full(n, 15.0)
    a_lumen_0 = np.full(n, 10.0)
    base_plaque_vol = np.clip(rng.normal(100.0, 15.0, size=n), 70.0, 140.0) # mm^3
    cur_plaque_rel = np.ones(n, dtype=np.float64) # 1.0 = 100%
    lipid_core_rel = np.ones(n, dtype=np.float64)
    t_cap = np.full(n, 54.0, dtype=np.float64)
    
    # Hemodynamics: baseline WSS in disturbed flow zone
    tau_w0 = np.clip(rng.normal(0.75, 0.15, size=n), 0.40, 1.10) # Pa
    cur_tau_w = tau_w0.copy()
    cur_ldl = base_ldl.copy()
    cur_lumen_area = a_lumen_0.copy()

    active = 0.0 if cfg.placebo_null_mode else 1.0
    inf_rate_per_step = (2.0 / 365.0) * cfg.dt

    for _ in range(steps):
        inf_spikes = rng.poisson(inf_rate_per_step, size=n) * 0.15

        # Bioactive concentrations
        ht_bio = active * adherence * (0.035 + p.k_deconj_macrophage * np.minimum(1.5, base_crp / 2.0))
        sac_bio = active * adherence * 1.8
        ua_bio = active * adherence * urolithin_factor * 1.2
        natto_bio = active * adherence * 0.4

        # eNOS activation
        raw_stim = (ht_bio * 2.5 + sac_bio * 0.8 + cur_tau_w * 1.0)
        enos_act = 1.0 + 0.70 * (raw_stim / (p.enos_hill_ec50 + raw_stim))

        # LDL & SREBP-2 feedback
        net_ldl_reduc = 0.14 * (1.0 - p.srebp2_rebound_gain) * active * adherence
        target_ldl = base_ldl * (1.0 - net_ldl_reduc)
        cur_ldl += (0.05 * cfg.dt) * (target_ldl - cur_ldl)

        # Lipid influx vs efflux
        oxldl_prot = 1.0 - 0.35 * np.minimum(1.0, ht_bio * 5.0)
        j_influx = 0.00030 * ((cur_ldl / 100.0) ** 1.5) * oxldl_prot
        rct_efflux = 0.00035 * (1.0 + 1.6 * (ua_bio / (0.8 + ua_bio)))

        d_lipid = (j_influx * 100.0 - rct_efflux * lipid_core_rel * 50.0) * cfg.dt
        lipid_core_rel = np.maximum(0.1, lipid_core_rel + (d_lipid / 50.0))

        # Macrophage MMP activity
        shear_sens = np.maximum(0.35, 1.15 - 0.5 * cur_tau_w)
        m1_drive = lipid_core_rel * shear_sens * (1.0 / enos_act)
        m2_drive = 0.4 * ua_bio + 0.3 * natto_bio + 0.2
        net_mmp = np.maximum(0.05, m1_drive - m2_drive) + (inf_spikes / 15.0)

        # GLAGOV & LUMEN COUPLING:
        # Plaque burden = (base_plaque * cur_plaque_rel) / A_EEM
        plaque_area_est = 5.0 * cur_plaque_rel + (t_cap - 54.0) * 0.02
        plaque_burden = plaque_area_est / a_eem_0
        
        # If burden < Glagov threshold: compensatory outward expansion preserves lumen
        # Beyond threshold: lumen narrows directly!
        lumen_encroachment = np.maximum(0.0, plaque_burden - p.glagov_threshold) * a_eem_0
        cur_lumen_area = np.maximum(1.5, a_lumen_0 - lumen_encroachment)
        
        # Poiseuille shear stress feedback: tau_w increases as lumen narrows!
        cur_tau_w = tau_w0 * ((a_lumen_0 / cur_lumen_area) ** 1.5)
        
        # Mechanosensitive feedback on VSMC collagen synthesis:
        # 1. High shear stress downregulates VSMC synthetic phenotype (anti-stenotic guard)
        wss_feedback = 1.0 / (1.0 + (cur_tau_w / p.wss_inhibition_ec50) ** 2)
        # 2. Diffusion barrier saturation: cap thickness cannot exceed hypoxia limit
        diffusion_guard = np.maximum(0.0, 1.0 - (t_cap / p.cap_max_diffusion))
        
        v_synth = p.vsmc_collagen_rate * (1.0 + 0.2 * active * adherence) * wss_feedback * diffusion_guard
        v_degrade = p.mmp_degradation_rate * net_mmp
        d_cap = (v_synth - v_degrade) * cfg.dt
        t_cap = np.maximum(15.0, t_cap + d_cap)

        # Update plaque volume
        d_plaque = (d_lipid * 0.85 + 0.003 * inf_spikes) * (cfg.dt / cfg.days * 730)
        cur_plaque_rel = np.maximum(0.1, cur_plaque_rel + (d_lipid * 0.85 / 100.0))

    final_plaque_abs_mm3 = base_plaque_vol * cur_plaque_rel
    delta_plaque_abs_mm3 = final_plaque_abs_mm3 - base_plaque_vol
    stenosis_pct = (1.0 - cur_lumen_area / a_lumen_0) * 100.0

    return {
        "plaque_pct_change": (cur_plaque_rel - 1.0) * 100.0,
        "delta_plaque_abs_mm3": delta_plaque_abs_mm3,
        "final_plaque_abs_mm3": final_plaque_abs_mm3,
        "t_cap_final": t_cap,
        "delta_cap_um": t_cap - 54.0,
        "stenosis_final_pct": stenosis_pct,
        "final_tau_w": cur_tau_w,
        "metabotype": metabotype,
    }

def run_coupled_evaluation():
    p = CoupledParams()
    cfg_super = CoupledConfig(n_patients=1000, seed=42)
    cfg_null = CoupledConfig(n_patients=1000, seed=42, placebo_null_mode=True)
    
    out_super = simulate_coupled_cohort(cfg_super, p)
    out_null = simulate_coupled_cohort(cfg_null, p)
    
    # Statistical synthesis
    def get_stats(arr):
        return {
            "median": round(float(np.median(arr)), 2),
            "ci_low": round(float(np.percentile(arr, 2.5)), 2),
            "ci_high": round(float(np.percentile(arr, 97.5)), 2),
            "iqr": round(float(np.percentile(arr, 75) - np.percentile(arr, 25)), 2)
        }
        
    results = {
        "Supermarket_Coupled": {
            "plaque_pct_change": get_stats(out_super["plaque_pct_change"]),
            "delta_plaque_abs_mm3": get_stats(out_super["delta_plaque_abs_mm3"]),
            "t_cap_final_um": get_stats(out_super["t_cap_final"]),
            "delta_cap_abs_um": get_stats(out_super["delta_cap_um"]),
            "stenosis_pct": get_stats(out_super["stenosis_final_pct"]),
            "final_tau_w_Pa": get_stats(out_super["final_tau_w"])
        },
        "Control_Null_Coupled": {
            "plaque_pct_change": get_stats(out_null["plaque_pct_change"]),
            "delta_plaque_abs_mm3": get_stats(out_null["delta_plaque_abs_mm3"]),
            "t_cap_final_um": get_stats(out_null["t_cap_final"]),
            "delta_cap_abs_um": get_stats(out_null["delta_cap_um"]),
            "stenosis_pct": get_stats(out_null["stenosis_final_pct"]),
            "final_tau_w_Pa": get_stats(out_null["final_tau_w"])
        }
    }
    
    # Print dual absolute & relative reporting
    print("\n=======================================================")
    print("COUPLED GLAGOV & LUMINAL SHEAR FEEDBACK EVALUATION")
    print("=======================================================")
    print("Dual Reporting: Absolute Shifts & Relative Metrics (N=1000, 24 Months):")
    print("\n1. PLAQUE-VOLUMEN (Absolut vs. Relativ):")
    print(f"  Kontrolle (Null):  Median {results['Control_Null_Coupled']['delta_plaque_abs_mm3']['median']:+.2f} mm3 [{results['Control_Null_Coupled']['delta_plaque_abs_mm3']['ci_low']:+.2f}, {results['Control_Null_Coupled']['delta_plaque_abs_mm3']['ci_high']:+.2f}]  | Relativ: {results['Control_Null_Coupled']['plaque_pct_change']['median']:+.2f}%")
    print(f"  Supermarkt:        Median {results['Supermarket_Coupled']['delta_plaque_abs_mm3']['median']:+.2f} mm3 [{results['Supermarket_Coupled']['delta_plaque_abs_mm3']['ci_low']:+.2f}, {results['Supermarket_Coupled']['delta_plaque_abs_mm3']['ci_high']:+.2f}]  | Relativ: {results['Supermarket_Coupled']['plaque_pct_change']['median']:+.2f}%")
    print(f"  -> Absolute Differenz: {results['Supermarket_Coupled']['delta_plaque_abs_mm3']['median'] - results['Control_Null_Coupled']['delta_plaque_abs_mm3']['median']:+.2f} mm3 Netto-Plaque-Vermeidung!")

    print("\n2. FIBRÖSE KAPPENDICKE & LUMEN-STENOSE (Glagov-Rückkopplung):")
    print(f"  Kontrolle (Null):  Enddicke: {results['Control_Null_Coupled']['t_cap_final_um']['median']:.1f} um (Delta: {results['Control_Null_Coupled']['delta_cap_abs_um']['median']:+.1f} um) | Stenose: {results['Control_Null_Coupled']['stenosis_pct']['median']:.1f}%")
    print(f"  Supermarkt:        Enddicke: {results['Supermarket_Coupled']['t_cap_final_um']['median']:.1f} um (Delta: {results['Supermarket_Coupled']['delta_cap_abs_um']['median']:+.1f} um) | Stenose: {results['Supermarket_Coupled']['stenosis_pct']['median']:.1f}%")
    print(f"  End-Scherspannung (tau_w): Kontrolle = {results['Control_Null_Coupled']['final_tau_w_Pa']['median']:.2f} Pa vs. Supermarkt = {results['Supermarket_Coupled']['final_tau_w_Pa']['median']:.2f} Pa")

    # ==============================================================================
    # RE-EVALUATE SOBOL SENSITIVITY UNDER LUMINAL COUPLING
    # ==============================================================================
    print("\nRunning Sobol re-evaluation under luminal feedback (SALib)...")
    from SALib.sample import sobol as sobol_sample
    from SALib.analyze import sobol as sobol_analyze

    problem = {
        "num_vars": 6,
        "names": [
            "k_deconj_macrophage",
            "srebp2_rebound_gain",
            "enos_hill_ec50",
            "vsmc_collagen_rate",
            "mmp_degradation_rate",
            "adherence_mean",
        ],
        "bounds": [
            [0.05, 0.70],
            [0.10, 0.35],
            [0.50, 2.50],
            [0.025, 0.065],
            [0.020, 0.055],
            [0.65, 0.90],
        ],
    }

    param_matrix = sobol_sample.sample(problem, N=128, calc_second_order=False) # 128*(6+2) = 1024 runs
    y_plaque_abs = np.zeros(len(param_matrix))
    y_cap_coupled = np.zeros(len(param_matrix))

    cfg_eval = CoupledConfig(n_patients=300, seed=123)
    for i, row in enumerate(param_matrix):
        hp = CoupledParams(
            k_deconj_macrophage=row[0],
            srebp2_rebound_gain=row[1],
            enos_hill_ec50=row[2],
            vsmc_collagen_rate=row[3],
            mmp_degradation_rate=row[4],
            adherence_mean=row[5]
        )
        out = simulate_coupled_cohort(cfg_eval, hp)
        y_plaque_abs[i] = np.median(out["delta_plaque_abs_mm3"])
        y_cap_coupled[i] = np.median(out["t_cap_final"])

    si_p = sobol_analyze.analyze(problem, y_plaque_abs, calc_second_order=False)
    si_c = sobol_analyze.analyze(problem, y_cap_coupled, calc_second_order=False)

    print("\nRe-Evaluated Sobol Indices under Luminal Coupling:")
    print(f"  {'Parameter':22s} | {'Plaque ST':10s} | {'Cap ST (Coupled)':16s}")
    print("  " + "-" * 55)
    for idx, name in enumerate(problem["names"]):
        st_p = float(si_p['ST'][idx])
        st_c = float(si_c['ST'][idx])
        print(f"  {name:22s} | {st_p:10.3f} | {st_c:16.3f}")

    results["sobol_coupled"] = {
        "names": problem["names"],
        "plaque_ST": [round(float(x), 4) for x in si_p['ST']],
        "cap_ST_coupled": [round(float(x), 4) for x in si_c['ST']]
    }

    # Plotting
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=300)
    
    # Left: Dual Absolute Plaque & Cap Shift Comparison
    ax1 = axs[0]
    cats = ['Kontrolle\n(Null)', 'Supermarkt\n(Gekoppelt)']
    delta_p_med = [results['Control_Null_Coupled']['delta_plaque_abs_mm3']['median'], results['Supermarket_Coupled']['delta_plaque_abs_mm3']['median']]
    delta_c_med = [results['Control_Null_Coupled']['delta_cap_abs_um']['median'], results['Supermarket_Coupled']['delta_cap_abs_um']['median']]
    
    x = np.arange(len(cats))
    w = 0.35
    ax1.bar(x - w/2, delta_p_med, w, label='Plaque-Änderung (mm³ absolut)', color='#d95f02', alpha=0.8)
    ax1.bar(x + w/2, delta_c_med, w, label='Kappen-Dickenänderung (µm absolut)', color='#1b9e77', alpha=0.8)
    for i in range(len(cats)):
        ax1.text(i - w/2, delta_p_med[i] + 0.5, f"{delta_p_med[i]:+.1f} mm³", ha='center', fontweight='bold', fontsize=9)
        ax1.text(i + w/2, delta_c_med[i] + 0.5, f"{delta_c_med[i]:+.1f} µm", ha='center', fontweight='bold', fontsize=9)
    ax1.set_xticks(x)
    ax1.set_xticklabels(cats, fontsize=11)
    ax1.set_ylabel("Absolute Änderung nach 2 Jahren", fontsize=11)
    ax1.set_title("A. Absolute Plaque- & Kappendynamik (Glagov-Lumen-Gegenkopplung)", fontsize=11, fontweight='bold')
    ax1.axhline(0, color='black', linestyle=':', lw=1.0)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9, loc='upper left')

    # Right: Re-evaluated Sobol Total-Order Indices
    ax2 = axs[1]
    clean_names = ["k_deconj", "srebp2", "enos_ec50", "vsmc_rate", "mmp_rate", "adherence"]
    st_c_vals = results["sobol_coupled"]["cap_ST_coupled"]
    bars = ax2.bar(clean_names, st_c_vals, color='#7570b3', alpha=0.85)
    for bar in bars:
        yval = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2.0, yval + 0.02, f"{yval:.2f}", ha='center', va='bottom', fontweight='bold', fontsize=9)
    ax2.set_title("B. Sobol ST auf Kappendicke (unter mechanosensitiver Koppelung)", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Totaler Sobol-Varianzanteil (ST)", fontsize=11)
    ax2.set_ylim(0, 1.15)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plot_path = os.path.join(RESULTS_DIR, "glagov_luminal_coupling_analysis.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Generated Glagov coupling plot: {plot_path}")

    with open(os.path.join(RESULTS_DIR, "glagov_luminal_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("Saved glagov_luminal_results.json successfully.")

if __name__ == "__main__":
    run_coupled_evaluation()
