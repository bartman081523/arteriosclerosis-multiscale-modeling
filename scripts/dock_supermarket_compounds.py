#!/usr/bin/env python3
"""
Download and dock supermarket bioactives against atheroprotective targets:
1. Human eNOS (1M9Q) - S-Allylcysteine & Hydroxytyrosol
2. Human Serum Paraoxonase-1 PON1 (1V04) - Urolithin A & Hydroxytyrosol
"""
import os
import subprocess
import json
import numpy as np
import requests
from Bio.PDB import PDBParser, PDBIO, Select
from vina import Vina

DATA_DIR = "/run/media/julian/ML6/arteriosklerosis/data"
RESULTS_DIR = "/run/media/julian/ML6/arteriosklerosis/results"
LIGAND_DIR = os.path.join(DATA_DIR, "ligands")

# 1. Download PDBs
urls = {
    "1M9Q.pdb": "https://files.rcsb.org/download/1M9Q.pdb", # Human eNOS
    "1V04.pdb": "https://files.rcsb.org/download/1V04.pdb"  # Serum Paraoxonase PON1
}

for fname, url in urls.items():
    fpath = os.path.join(DATA_DIR, fname)
    if not os.path.exists(fpath):
        r = requests.get(url, timeout=30)
        with open(fpath, "wb") as f:
            f.write(r.content)
        print(f"Downloaded {fname}")

# 2. Clean and convert to PDBQT
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == ' '

parser = PDBParser(QUIET=True)
io = PDBIO()

# eNOS
s_enos = parser.get_structure("eNOS", os.path.join(DATA_DIR, "1M9Q.pdb"))
clean_enos_pdb = os.path.join(DATA_DIR, "enos_clean.pdb")
io.set_structure(s_enos)
io.save(clean_enos_pdb, ProteinSelect())
enos_pdbqt = os.path.join(DATA_DIR, "enos_clean.pdbqt")
subprocess.run(["obabel", "-ipdb", clean_enos_pdb, "-opdbqt", "-O", enos_pdbqt, "-xr", "-xh"], check=True)

# PON1
s_pon1 = parser.get_structure("PON1", os.path.join(DATA_DIR, "1V04.pdb"))
clean_pon1_pdb = os.path.join(DATA_DIR, "pon1_clean.pdb")
io.set_structure(s_pon1)
io.save(clean_pon1_pdb, ProteinSelect())
pon1_pdbqt = os.path.join(DATA_DIR, "pon1_clean.pdbqt")
subprocess.run(["obabel", "-ipdb", clean_pon1_pdb, "-opdbqt", "-O", pon1_pdbqt, "-xr", "-xh"], check=True)

# Sanitize PDBQTs
for p in [enos_pdbqt, pon1_pdbqt]:
    with open(p) as f:
        lines = f.readlines()
    clean_lines = [l for l in lines if l.startswith(("ATOM", "HETATM", "TER", "END"))]
    with open(p, "w") as f:
        f.writelines(clean_lines)

# 3. Calculate grid box centers
def get_centroid(pdb_path):
    coords = []
    with open(pdb_path) as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.mean(np.array(coords), axis=0)

enos_center = get_centroid(enos_pdbqt)
pon1_center = get_centroid(pon1_pdbqt)

print(f"eNOS Center: {enos_center}")
print(f"PON1 Center: {pon1_center}")

# 4. Docking runs
dockings = [
    ("eNOS", enos_pdbqt, enos_center, "s_allylcysteine", "Knoblauch / S-Allylcystein an eNOS"),
    ("eNOS", enos_pdbqt, enos_center, "hydroxytyrosol", "Olivenöl / Hydroxytyrosol an eNOS"),
    ("PON1", pon1_pdbqt, pon1_center, "urolithin_a", "Granatapfel / Urolithin A an PON1"),
    ("PON1", pon1_pdbqt, pon1_center, "hydroxytyrosol", "Olivenöl / Hydroxytyrosol an PON1")
]

results = {}
for target_name, rec_path, center, lig_name, desc in dockings:
    v = Vina(sf_name="vina", cpu=4, seed=42)
    v.set_receptor(rec_path)
    lig_path = os.path.join(LIGAND_DIR, f"{lig_name}.pdbqt")
    v.set_ligand_from_file(lig_path)
    v.compute_vina_maps(center=center.tolist(), box_size=[22.0, 22.0, 22.0])
    v.dock(exhaustiveness=16, n_poses=3)
    energies = v.energies(n_poses=3)
    dG = float(energies[0][0])
    kd_uM = np.exp(dG / (1.9872e-3 * 310.15)) * 1e6
    out_pose = os.path.join(RESULTS_DIR, f"{target_name}_{lig_name}_docked.pdbqt")
    v.write_poses(out_pose, n_poses=3, overwrite=True)
    
    key = f"{target_name}_{lig_name}"
    results[key] = {
        "target": target_name,
        "ligand": lig_name,
        "description": desc,
        "affinity_kcal_mol": round(dG, 2),
        "est_kd_uM": round(float(kd_uM), 2)
    }
    print(f"{desc:45} | Affinity: {dG:6.2f} kcal/mol | Est. Kd: {kd_uM:8.2f} uM")

with open(os.path.join(RESULTS_DIR, "supermarket_docking_results.json"), "w") as f:
    json.dump(results, f, indent=2)

print("\nSaved supermarket_docking_results.json successfully.")
