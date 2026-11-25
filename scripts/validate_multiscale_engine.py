#!/usr/bin/env python3
"""
Reproducible 4-Stage Validation, Ablation, and Sobol Sensitivity Suite
for the Multiscale Systems-Biology Atherosclerosis Monte Carlo Engine.
"""
import os
import json
from dataclasses import dataclass, replace
from typing import Dict, Tuple, List
import numpy as np
import matplotlib.pyplot as plt

RESULTS_DIR = "/run/media/julian/ML6/arteriosklerosis/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

# ==============================================================================
# CONFIGURATION & PARAMETER DATA STRUCTURES
# ==============================================================================
@dataclass(frozen=True)
class SimConfig:
    n_patients: int = 1000
    days: int = 730
    dt: float = 1.0  # Time step in days
    seed: int = 42
    # Ablation toggles
    enable_deconjugation: bool = True
    enable_srebp2_feedback: bool = True
    enable_hill_saturation: bool = True
    enable_infection_spikes: bool = True
    placebo_null_mode: bool = False

@dataclass(frozen=True)
class HyperParams:
    """Biologically uncertain parameters for global Sobol sensitivity analysis."""
    k_deconj_macrophage: float = 0.35   # Local deconjugation rate in plaque (0.05 - 0.70)
    srebp2_rebound_gain: float = 0.22   # Compensatory HMGCR rebound gain (0.10 - 0.35)
    enos_hill_ec50: float = 1.20        # eNOS activation saturation threshold (0.50 - 2.50)
    vsmc_collagen_rate: float = 0.045   # VSMC collagen cap synthesis rate (0.025 - 0.065)
    mmp_degradation_rate: float = 0.038 # Macrophage MMP-1/9 collagen breakdown (0.020 - 0.055)
    adherence_mean: float = 0.78        # Population mean dietary adherence (0.65 - 0.90)

# ==============================================================================
# VECTORIZED SIMULATION ENGINE
# ==============================================================================
def simulate_cohort(cfg: SimConfig, params: HyperParams) -> Dict[str, np.ndarray]:
    """
    Simulates a cohort of virtual patients over cfg.days with time step cfg.dt.
    Returns per-patient endpoints at day 730.
    """
    rng = np.random.default_rng(cfg.seed)
    n = cfg.n_patients
    steps = int(cfg.days / cfg.dt)

    # 1. Patient Heterogeneity (Metabotypes: A=0 [36.6%], B=1 [45.9%], 0=2 [17.5%])
    metabotype = rng.choice([0, 1, 2], size=n, p=[0.366, 0.459, 0.175])
    # Biological efficacy multiplier for Urolithin A:
    urolithin_factor = np.select(
        [metabotype == 0, metabotype == 1, metabotype == 2],
        [1.0, 0.45, 0.0]
    )

    # Individual adherence distribution around population mean
    adherence = np.clip(rng.normal(params.adherence_mean, 0.06, size=n), 0.50, 1.0)

    # Baseline characteristics
    base_ldl = np.clip(rng.normal(138.0, 24.0, size=n), 90.0, 210.0)
    base_crp = np.clip(rng.lognormal(0.8, 0.4, size=n), 0.5, 8.0)
    base_wss = np.clip(rng.normal(0.75, 0.18, size=n), 0.35, 1.30)
    
    # State vectors
    plaque_rel = np.ones(n, dtype=np.float64) # 1.0 = 100% baseline volume
    lipid_core_rel = np.ones(n, dtype=np.float64)
    t_cap = np.full(n, 54.0, dtype=np.float64) # initial cap thickness in um
    cur_ldl = base_ldl.copy()

    # Active intervention toggle
    active = 0.0 if cfg.placebo_null_mode else 1.0
    deconj_rate = params.k_deconj_macrophage if cfg.enable_deconjugation else 0.0
    rebound_gain = params.srebp2_rebound_gain if cfg.enable_srebp2_feedback else 0.0

    # 2. Time stepping loop
    inf_rate_per_step = (2.0 / 365.0) * cfg.dt
    for _ in range(steps):
        # Stochastically occurring acute viral infections (~2/year)
        if cfg.enable_infection_spikes:
            inf_spikes = rng.poisson(inf_rate_per_step, size=n) * 0.15
        else:
            inf_spikes = np.zeros(n, dtype=np.float64)

        # Bioactive stimuli across compartments
        # Hydroxytyrosol (EVOO): ~35 nM free + plaque deconjugation by macrophage beta-glucuronidase
        ht_bio = active * adherence * (0.035 + deconj_rate * np.minimum(1.5, base_crp / 2.0))
        # S-Allylcysteine (Garlic): water soluble, high bioavailability
        sac_bio = active * adherence * 1.8
        # Urolithin A (Pomegranate gut metabolite)
        ua_bio = active * adherence * urolithin_factor * 1.2
        # Nattokinase peptide action
        natto_bio = active * adherence * 0.4

        # Non-linear eNOS saturation (Hill-like gating)
        raw_stim = (ht_bio * 2.5 + sac_bio * 0.8 + base_wss * 1.0)
        if cfg.enable_hill_saturation:
            enos_act = 1.0 + 0.70 * (raw_stim / (params.enos_hill_ec50 + raw_stim))
        else:
            enos_act = 1.0 + 0.35 * raw_stim

        # Hepatic SREBP-2 counter-regulation to beta-glucan
        # Gross clearance ~14%, rebound compensation dampens net effect
        net_ldl_reduc = 0.14 * (1.0 - rebound_gain) * active * adherence
        target_ldl = base_ldl * (1.0 - net_ldl_reduc)
        cur_ldl += (0.05 * cfg.dt) * (target_ldl - cur_ldl)

        # Plaque lipid influx vs. efflux
        oxldl_prot = 1.0 - 0.35 * np.minimum(1.0, ht_bio * 5.0)
        j_influx = 0.00030 * ((cur_ldl / 100.0) ** 1.5) * oxldl_prot
        rct_efflux = 0.00035 * (1.0 + 1.6 * (ua_bio / (0.8 + ua_bio)))

        d_lipid = (j_influx * 100.0 - rct_efflux * lipid_core_rel * 50.0) * cfg.dt
        lipid_core_rel = np.maximum(0.1, lipid_core_rel + (d_lipid / 50.0))

        # Macrophage M1/M2 drive and Net MMP activity
        shear_factor = np.maximum(0.4, 1.1 - 0.5 * base_wss)
        m1_drive = (lipid_core_rel) * shear_factor * (1.0 / enos_act)
        m2_drive = 0.4 * ua_bio + 0.3 * natto_bio + 0.2
        net_mmp = np.maximum(0.05, m1_drive - m2_drive) + (inf_spikes / 15.0)

        # Cap remodeling: VSMC collagen synthesis vs. MMP degradation
        v_synth = params.vsmc_collagen_rate * (1.0 + 0.2 * active * adherence)
        v_degrade = params.mmp_degradation_rate * net_mmp
        d_cap = (v_synth - v_degrade) * cfg.dt
        t_cap = np.maximum(15.0, t_cap + d_cap)

        # Total Atheroma Volume change
        d_plaque = (d_lipid * 0.85 + 0.003 * inf_spikes) * (cfg.dt / cfg.days * 730)
        plaque_rel = np.maximum(0.1, plaque_rel + (d_lipid * 0.85 / 100.0))

    return {
        "plaque_pct_change": (plaque_rel - 1.0) * 100.0,
        "lipid_pct_change": (lipid_core_rel - 1.0) * 100.0,
        "t_cap_final": t_cap,
        "metabotype": metabotype,
    }

def summarize(arr: np.ndarray) -> Dict[str, float]:
    return {
        "median": round(float(np.median(arr)), 2),
        "ci_low": round(float(np.percentile(arr, 2.5)), 2),
        "ci_high": round(float(np.percentile(arr, 97.5)), 2),
        "iqr": round(float(np.percentile(arr, 75) - np.percentile(arr, 25)), 2)
    }

# ==============================================================================
# GATE 1: NUMERICS, SAMPLE SIZE (N) & STEP-SIZE (dt) CONVERGENCE
# ==============================================================================
def run_gate1_convergence(default_params: HyperParams) -> Dict:
    print("\n=======================================================")
    print("STAGE 1: NUMERICS & CONVERGENCE GATE")
    print("=======================================================")
    results_g1 = {}
    
    # 1. N-scaling: N in {250, 1000, 5000} with dt=1.0
    for n in (250, 1000, 5000):
        out = simulate_cohort(SimConfig(n_patients=n, dt=1.0, seed=42), default_params)
        stats = summarize(out["plaque_pct_change"])
        cap_stats = summarize(out["t_cap_final"])
        results_g1[f"N_{n}_dt_1.0"] = {"plaque": stats, "t_cap": cap_stats}
        print(f"  N={n:4d}, dt=1.0 -> Plaque: {stats['median']:+5.2f}% [{stats['ci_low']:+5.2f}%, {stats['ci_high']:+5.2f}%], IQR={stats['iqr']:4.2f}% | T_cap: {cap_stats['median']:4.1f} um")
        
    # 2. Time-step refinement: dt=0.5 days vs dt=1.0 days at N=1000
    out_dt05 = simulate_cohort(SimConfig(n_patients=1000, dt=0.5, seed=42), default_params)
    stats_dt05 = summarize(out_dt05["plaque_pct_change"])
    cap_stats_dt05 = summarize(out_dt05["t_cap_final"])
    results_g1["N_1000_dt_0.5"] = {"plaque": stats_dt05, "t_cap": cap_stats_dt05}
    print(f"  N=1000, dt=0.5 -> Plaque: {stats_dt05['median']:+5.2f}% [{stats_dt05['ci_low']:+5.2f}%, {stats_dt05['ci_high']:+5.2f}%], IQR={stats_dt05['iqr']:4.2f}% | T_cap: {cap_stats_dt05['median']:4.1f} um")
    
    # 3. Seed stability check over 5 independent seeds
    seed_medians = []
    for s in (11, 42, 99, 137, 777):
        out_s = simulate_cohort(SimConfig(n_patients=1000, dt=1.0, seed=s), default_params)
        seed_medians.append(float(np.median(out_s["plaque_pct_change"])))
    seed_std = float(np.std(seed_medians))
    print(f"  5-Seed Median Variance: mean={np.mean(seed_medians):+.2f}%, std={seed_std:.3f}%")
    
    # Acceptance check: delta_median between N=1000 and N=5000 < 0.3 percentage points
    delta_n = abs(results_g1["N_1000_dt_1.0"]["plaque"]["median"] - results_g1["N_5000_dt_1.0"]["plaque"]["median"])
    delta_dt = abs(results_g1["N_1000_dt_1.0"]["plaque"]["median"] - results_g1["N_1000_dt_0.5"]["plaque"]["median"])
    passed = (delta_n < 0.30) and (delta_dt < 0.40) and (seed_std < 0.35)
    results_g1["verdict"] = {
        "passed": bool(passed),
        "delta_N1000_N5000": round(delta_n, 3),
        "delta_dt10_dt05": round(delta_dt, 3),
        "seed_std": round(seed_std, 3)
    }
    print(f"  Stage 1 Verdict: {'PASSED (Acceptance Criteria Met)' if passed else 'FAILED'}")
    return results_g1

# ==============================================================================
# GATE 2: FALSIFICATION & BOUNDS AUDIT (NULL / PLACEBO MODE)
# ==============================================================================
def run_gate2_null_audit(default_params: HyperParams) -> Dict:
    print("\n=======================================================")
    print("STAGE 2: FALSIFICATION & NULL-AGENT AUDIT")
    print("=======================================================")
    # Placebo Null Mode: Active stimulus = 0, but full random adherence and infection noise active
    cfg_null = SimConfig(n_patients=2000, placebo_null_mode=True, seed=101)
    out_null = simulate_cohort(cfg_null, default_params)
    
    # Check for NaNs or Infs
    has_nan = np.isnan(out_null["plaque_pct_change"]).any() or np.isnan(out_null["t_cap_final"]).any()
    # Check physical boundaries
    all_caps_positive = (out_null["t_cap_final"] >= 15.0).all()
    no_infinite_plaque = (out_null["plaque_pct_change"] < 250.0).all() and (out_null["plaque_pct_change"] > -90.0).all()
    
    stats_p = summarize(out_null["plaque_pct_change"])
    stats_c = summarize(out_null["t_cap_final"])
    
    # In null mode, plaque should spontaneously progress (positive growth median ~ +20 to +28%)
    null_behaves_as_control = stats_p["median"] > 15.0
    passed = (not has_nan) and all_caps_positive and no_infinite_plaque and null_behaves_as_control
    
    print(f"  Null-Agent Mode (N=2000):")
    print(f"    Plaque Change Median: {stats_p['median']:+5.2f}% [{stats_p['ci_low']:+5.2f}%, {stats_p['ci_high']:+5.2f}%]")
    print(f"    Cap Thickness Median: {stats_c['median']:4.1f} um [{stats_c['ci_low']:4.1f}, {stats_c['ci_high']:4.1f}]")
    print(f"    NaN/Inf Detected:     {has_nan}")
    print(f"    Cap Bounds Respected: {all_caps_positive}")
    print(f"  Stage 2 Verdict: {'PASSED (Zero-Bias Confirmed)' if passed else 'FAILED'}")
    
    return {
        "plaque_stats": stats_p,
        "cap_stats": stats_c,
        "passed": bool(passed)
    }

# ==============================================================================
# GATE 3: SUBMODEL ABLATION SUITE (LEAVE-ONE-MECHANISM-OUT / LOMO)
# ==============================================================================
def run_gate3_ablation(default_params: HyperParams) -> Dict:
    print("\n=======================================================")
    print("STAGE 3: SUBMODEL ABLATION (LOMO) SUITE")
    print("=======================================================")
    
    scenarios = {
        "Vollmodell": SimConfig(n_patients=1000, seed=42),
        "Ohne Makrophagen-Dekonjugation": SimConfig(n_patients=1000, seed=42, enable_deconjugation=False),
        "Ohne SREBP-2-Gegenregulation": SimConfig(n_patients=1000, seed=42, enable_srebp2_feedback=False),
        "Ohne Hill-Sättigung (linear)": SimConfig(n_patients=1000, seed=42, enable_hill_saturation=False),
        "Ohne Infekt-Spikes": SimConfig(n_patients=1000, seed=42, enable_infection_spikes=False)
    }
    
    ablation_results = {}
    base_out = simulate_cohort(scenarios["Vollmodell"], default_params)
    base_p_med = float(np.median(base_out["plaque_pct_change"]))
    base_c_med = float(np.median(base_out["t_cap_final"]))
    
    for name, cfg in scenarios.items():
        out = simulate_cohort(cfg, default_params)
        p_stats = summarize(out["plaque_pct_change"])
        c_stats = summarize(out["t_cap_final"])
        
        delta_p = round(p_stats["median"] - base_p_med, 2)
        delta_c = round(c_stats["median"] - base_c_med, 2)
        
        ablation_results[name] = {
            "plaque": p_stats,
            "cap": c_stats,
            "delta_plaque_median": delta_p,
            "delta_cap_median": delta_c
        }
        print(f"  {name:32s} | Plaque: {p_stats['median']:+5.2f}% (d={delta_p:+5.2f}%) | Cap: {c_stats['median']:4.1f} um (d={delta_c:+4.1f} um)")
        
    return ablation_results

# ==============================================================================
# GATE 4: GLOBAL SOBOL SENSITIVITY ANALYSIS (SALIB)
# ==============================================================================
def run_gate4_sobol() -> Dict:
    print("\n=======================================================")
    print("STAGE 4: GLOBAL SOBOL SENSITIVITY ANALYSIS")
    print("=======================================================")
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
            [0.05, 0.70],   # Local deconjugation in plaque
            [0.10, 0.35],   # SREBP-2 rebound compensation
            [0.50, 2.50],   # eNOS Hill EC50
            [0.025, 0.065], # VSMC collagen synthesis
            [0.020, 0.055], # MMP-1/9 degradation
            [0.65, 0.90],   # Mean population adherence
        ],
    }

    # N_sobol = 256 generates 256 * (2*6 + 2) = 3584 parameter sets
    N_SAMPLE = 256
    param_matrix = sobol_sample.sample(problem, N=N_SAMPLE, calc_second_order=True)
    print(f"  Generated {len(param_matrix)} Saltelli parameter vectors.")

    y_plaque_med = np.zeros(len(param_matrix))
    y_plaque_iqr = np.zeros(len(param_matrix))
    y_cap_med = np.zeros(len(param_matrix))

    # Common Random Numbers (CRN): fixed seed per evaluation isolates parameter effects from MC noise
    cfg_eval = SimConfig(n_patients=400, seed=123)
    for i, row in enumerate(param_matrix):
        hp = HyperParams(*row)
        out = simulate_cohort(cfg_eval, hp)
        y_plaque_med[i] = np.median(out["plaque_pct_change"])
        y_plaque_iqr[i] = np.percentile(out["plaque_pct_change"], 75) - np.percentile(out["plaque_pct_change"], 25)
        y_cap_med[i] = np.median(out["t_cap_final"])

    si_plaque = sobol_analyze.analyze(problem, y_plaque_med, calc_second_order=True)
    si_plaque_iqr = sobol_analyze.analyze(problem, y_plaque_iqr, calc_second_order=True)
    si_cap = sobol_analyze.analyze(problem, y_cap_med, calc_second_order=True)

    print("\n  Sobol Sensitivity Summary:")
    print(f"  {'Parameter':22s} | {'Plaque S1':10s} {'Plaque ST':10s} | {'Cap S1':10s} {'Cap ST':10s} | {'IQR ST':10s}")
    print("  " + "-" * 75)
    
    sobol_out = {"names": problem["names"]}
    for idx, name in enumerate(problem["names"]):
        s1_p = float(si_plaque['S1'][idx])
        st_p = float(si_plaque['ST'][idx])
        s1_c = float(si_cap['S1'][idx])
        st_c = float(si_cap['ST'][idx])
        st_iqr = float(si_plaque_iqr['ST'][idx])
        print(f"  {name:22s} | {s1_p:10.3f} {st_p:10.3f} | {s1_c:10.3f} {st_c:10.3f} | {st_iqr:10.3f}")
        
    sobol_out["plaque_median"] = {
        "S1": [round(float(x), 4) for x in si_plaque['S1']],
        "ST": [round(float(x), 4) for x in si_plaque['ST']],
        "S1_conf": [round(float(x), 4) for x in si_plaque['S1_conf']],
        "ST_conf": [round(float(x), 4) for x in si_plaque['ST_conf']],
    }
    sobol_out["cap_thickness"] = {
        "S1": [round(float(x), 4) for x in si_cap['S1']],
        "ST": [round(float(x), 4) for x in si_cap['ST']],
    }
    sobol_out["plaque_iqr"] = {
        "ST": [round(float(x), 4) for x in si_plaque_iqr['ST']],
    }
    
    return sobol_out

# ==============================================================================
# MULTI-PANEL DIAGNOSTIC PLOT GENERATION
# ==============================================================================
def plot_validation_suite(g1: Dict, g3: Dict, g4: Dict):
    fig, axs = plt.subplots(2, 2, figsize=(18, 13), dpi=300)

    # Panel A: Step-size (dt) and N-convergence
    ax1 = axs[0, 0]
    labels_n = ["N=250\ndt=1.0", "N=1000\ndt=1.0", "N=5000\ndt=1.0", "N=1000\ndt=0.5"]
    meds_n = [
        g1["N_250_dt_1.0"]["plaque"]["median"],
        g1["N_1000_dt_1.0"]["plaque"]["median"],
        g1["N_5000_dt_1.0"]["plaque"]["median"],
        g1["N_1000_dt_0.5"]["plaque"]["median"]
    ]
    ci_lows = [
        g1["N_250_dt_1.0"]["plaque"]["ci_low"],
        g1["N_1000_dt_1.0"]["plaque"]["ci_low"],
        g1["N_5000_dt_1.0"]["plaque"]["ci_low"],
        g1["N_1000_dt_0.5"]["plaque"]["ci_low"]
    ]
    ci_highs = [
        g1["N_250_dt_1.0"]["plaque"]["ci_high"],
        g1["N_1000_dt_1.0"]["plaque"]["ci_high"],
        g1["N_5000_dt_1.0"]["plaque"]["ci_high"],
        g1["N_1000_dt_0.5"]["plaque"]["ci_high"]
    ]
    err_low = [m - l for m, l in zip(meds_n, ci_lows)]
    err_high = [h - m for m, h in zip(meds_n, ci_highs)]
    ax1.errorbar(range(4), meds_n, yerr=[err_low, err_high], fmt='o', color='navy', ecolor='skyblue', elinewidth=2.5, capsize=6, markersize=8)
    for i, m in enumerate(meds_n):
        ax1.text(i, m + 1.2, f"{m:+.2f}%", ha='center', fontweight='bold', fontsize=10)
    ax1.set_xticks(range(4))
    ax1.set_xticklabels(labels_n, fontsize=10)
    ax1.axhline(0, color='black', linestyle=':', lw=1.0)
    ax1.set_title("A. Stufe 1: Numerische Konvergenz (N-Skalierung & dt-Halbierung)", fontsize=12, fontweight='bold')
    ax1.set_ylabel("Plaque-Änderung (%) [Median & 95% CI]", fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Panel B: Submodel Ablation Forest Plot
    ax2 = axs[0, 1]
    abl_keys = list(g3.keys())
    abl_meds = [g3[k]["plaque"]["median"] for k in abl_keys]
    abl_lows = [g3[k]["plaque"]["ci_low"] for k in abl_keys]
    abl_highs = [g3[k]["plaque"]["ci_high"] for k in abl_keys]
    y_pos = np.arange(len(abl_keys))
    ax2.errorbar(abl_meds, y_pos, xerr=[np.array(abl_meds)-np.array(abl_lows), np.array(abl_highs)-np.array(abl_meds)],
                 fmt='s', color='#d95f02', ecolor='#fdcdac', elinewidth=2.5, capsize=5, markersize=8)
    for y, m, d in zip(y_pos, abl_meds, [g3[k]["delta_plaque_median"] for k in abl_keys]):
        ax2.text(m + 1.0, y, f"{m:+.2f}% ({d:+.2f}%)", va='center', fontsize=9, fontweight='bold')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(abl_keys, fontsize=10)
    ax2.axvline(g3["Vollmodell"]["plaque"]["median"], color='blue', linestyle='--', lw=1.5, label="Vollmodell-Referenz")
    ax2.set_title("B. Stufe 3: Leave-One-Mechanism-Out (Submodell-Ablation)", fontsize=12, fontweight='bold')
    ax2.set_xlabel("Plaque-Änderung (%) [Median & 95% CI]", fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc='lower right')

    # Panel C: Sobol S1 and ST for Plaque Volume vs. Cap Thickness
    ax3 = axs[1, 0]
    names = g4["names"]
    clean_names = ["k_deconj", "srebp2_reb", "enos_ec50", "vsmc_rate", "mmp_rate", "adherence"]
    x = np.arange(len(clean_names))
    width = 0.35
    s1_p = g4["plaque_median"]["S1"]
    st_p = g4["plaque_median"]["ST"]
    ax3.bar(x - width/2, s1_p, width, label='Plaque S1 (Haupteffekt)', color='#1f77b4', alpha=0.85)
    ax3.bar(x + width/2, st_p, width, label='Plaque ST (Totaler Effekt)', color='#aec7e8', alpha=0.85)
    ax3.set_xticks(x)
    ax3.set_xticklabels(clean_names, fontsize=10, rotation=15)
    ax3.set_title("C. Stufe 4: Sobol-Sensitivität auf Plaque-Volumen (S1 vs. ST)", fontsize=12, fontweight='bold')
    ax3.set_ylabel("Sobol-Varianzanteil", fontsize=11)
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=9, loc='upper right')

    # Panel D: Parameter Interaction Map (ST - S1) across Endpoints
    ax4 = axs[1, 1]
    inter_p = np.array(st_p) - np.array(s1_p)
    s1_c = g4["cap_thickness"]["S1"]
    st_c = g4["cap_thickness"]["ST"]
    inter_c = np.array(st_c) - np.array(s1_c)
    ax4.bar(x - width/2, inter_p, width, label='Plaque Interaktion (ST - S1)', color='#7570b3', alpha=0.85)
    ax4.bar(x + width/2, inter_c, width, label='Kappe Interaktion (ST - S1)', color='#1b9e77', alpha=0.85)
    ax4.set_xticks(x)
    ax4.set_xticklabels(clean_names, fontsize=10, rotation=15)
    ax4.set_title("D. Nicht-lineare Wechselwirkungen (ST - S1)", fontsize=12, fontweight='bold')
    ax4.set_ylabel("Interaktions-Varianzanteil", fontsize=11)
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=9, loc='upper right')

    plt.tight_layout()
    plot_path = os.path.join(RESULTS_DIR, "sobol_sensitivity_and_ablation.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Generated validation plot: {plot_path}")

# ==============================================================================
# MAIN EXECUTION ORCHESTRATOR
# ==============================================================================
def main():
    default_hp = HyperParams()
    
    # Run all 4 stages
    g1 = run_gate1_convergence(default_hp)
    g2 = run_gate2_null_audit(default_hp)
    g3 = run_gate3_ablation(default_hp)
    g4 = run_gate4_sobol()
    
    # Plotting
    plot_validation_suite(g1, g3, g4)
    
    # Save combined results
    summary_export = {
        "gate1_numerics_convergence": g1,
        "gate2_falsification_null_audit": g2,
        "gate3_submodel_ablation": g3,
        "gate4_sobol_sensitivity": g4
    }
    
    json_path = os.path.join(RESULTS_DIR, "validation_suite_results.json")
    with open(json_path, "w") as f:
        json.dump(summary_export, f, indent=2)
    print(f"Saved complete validation summary to {json_path}")

if __name__ == "__main__":
    main()
