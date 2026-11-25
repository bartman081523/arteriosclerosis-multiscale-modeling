#!/usr/bin/env python3
"""
Molecular Docking of Acetate and Short-Chain Fatty Acids against:
1. Human Free Fatty Acid Receptor 2 (FFAR2 / GPR43, PDB 8J24)
2. Human Acetyl-CoA Synthetase 2 (ACSS2, AlphaFold AF-Q9NR19-F1)
"""
import os
import json
import math
import numpy as np
from vina import Vina
from Bio.PDB import PDBParser

BASE_DIR = "/run/media/julian/ML6/arteriosklerosis"
DATA_DIR = os.path.join(BASE_DIR, "data")
LIGAND_DIR = os.path.join(DATA_DIR, "ligands")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

with open(os.path.join(DATA_DIR, "grid_boxes.json")) as f:
    grid_boxes = json.load(f)

targets = {
    "FFAR2": {
        "receptor_pdbqt": os.path.join(DATA_DIR, "ffar2_clean.pdbqt"),
        "box": grid_boxes["ffar2"],
        "description": "Human Free Fatty Acid Receptor 2 (GPR43, Cryo-EM 8J24, 2.6A)"
    },
    "ACSS2": {
        "receptor_pdbqt": os.path.join(DATA_DIR, "acss2_clean.pdbqt"),
        "box": grid_boxes["acss2"],
        "description": "Human Acetyl-CoA Synthetase 2 (AlphaFold AF-Q9NR19-F1)"
    }
}

ligands = ["acetate", "acetic_acid_neutral", "formate", "propionate", "butyrate"]

# Gas constant R = 1.9872 cal/(mol*K), T = 310.15 K (37 C)
R = 1.9872e-3  # kcal/(mol*K)
T = 310.15     # Kelvin

docking_results = {}

for tname, tinfo in targets.items():
    print(f"\n=======================================================")
    print(f"Starting Docking for Target: {tname}")
    print(f"Description: {tinfo['description']}")
    print(f"=======================================================")
    
    docking_results[tname] = {}
    box = tinfo["box"]
    
    for lig in ligands:
        lig_path = os.path.join(LIGAND_DIR, f"{lig}.pdbqt")
        out_path = os.path.join(RESULTS_DIR, f"{tname}_{lig}_docked.pdbqt")
        
        v = Vina(sf_name="vina", cpu=4, seed=42)
        v.set_receptor(tinfo["receptor_pdbqt"])
        v.set_ligand_from_file(lig_path)
        v.compute_vina_maps(
            center=[box["center_x"], box["center_y"], box["center_z"]],
            box_size=[box["size_x"], box["size_y"], box["size_z"]]
        )
        
        # Dock with exhaustiveness = 16 for robust sampling of small carboxylate ligands
        v.dock(exhaustiveness=16, n_poses=5)
        energies = v.energies(n_poses=5)
        best_affinity = float(energies[0][0])  # kcal/mol
        
        # Estimate Kd (M) from deltaG = R*T*ln(Kd)  => Kd = exp(deltaG / (R*T))
        kd_molar = math.exp(best_affinity / (R * T))
        kd_micromolar = kd_molar * 1e6
        
        v.write_poses(out_path, n_poses=5, overwrite=True)
        
        docking_results[tname][lig] = {
            "affinity_kcal_mol": round(best_affinity, 2),
            "estimated_kd_uM": round(kd_micromolar, 2),
            "all_energies": [round(float(e[0]), 2) for e in energies],
            "docked_pdbqt": out_path
        }
        print(f"  Ligand: {lig:20} | Affinity: {best_affinity:6.2f} kcal/mol | Est. Kd: {kd_micromolar:10.2f} uM")

# Save results to JSON
res_json_path = os.path.join(RESULTS_DIR, "docking_results.json")
with open(res_json_path, "w") as f:
    json.dump(docking_results, f, indent=2)

print(f"\nSaved docking results to {res_json_path}")
