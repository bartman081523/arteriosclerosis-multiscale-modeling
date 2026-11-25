#!/usr/bin/env python3
"""
Extract receptor chains and convert to PDBQT for docking.
Calculate docking grid box centers for FFAR2 (8J24) and ACSS2 (AF-Q9NR19).
"""
import os
import subprocess
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

DATA_DIR = "/run/media/julian/ML6/arteriosklerosis/data"

class FFAR2Select(Select):
    def accept_chain(self, chain):
        return chain.id == 'D'
    def accept_residue(self, residue):
        # Exclude heteroatoms / ligands / waters
        return residue.id[0] == ' '

# 1. FFAR2 extraction
parser = PDBParser(QUIET=True)
s_8j24 = parser.get_structure("8J24", os.path.join(DATA_DIR, "8J24.pdb"))
io = PDBIO()
io.set_structure(s_8j24)
ffar2_clean_pdb = os.path.join(DATA_DIR, "ffar2_clean.pdb")
io.save(ffar2_clean_pdb, FFAR2Select())
print(f"Saved clean FFAR2 to {ffar2_clean_pdb}")

# Convert FFAR2 to PDBQT
ffar2_pdbqt = os.path.join(DATA_DIR, "ffar2_clean.pdbqt")
subprocess.run([
    "obabel", "-ipdb", ffar2_clean_pdb,
    "-opdbqt", "-O", ffar2_pdbqt,
    "-xr", "-xh"
], check=True)
print(f"Generated {ffar2_pdbqt}")

# Grid center for FFAR2 from native ACT ligand coords:
# C: 108.150, 114.982, 128.566
ffar2_box = {
    "center_x": 108.15,
    "center_y": 114.98,
    "center_z": 128.57,
    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0
}
print("FFAR2 Docking Grid Box:", ffar2_box)

# 2. ACSS2 preparation
s_acss2 = parser.get_structure("ACSS2", os.path.join(DATA_DIR, "AF-Q9NR19-F1.pdb"))
# Active site residues: Arg202, Thr315, Lys363
active_coords = []
for model in s_acss2:
    for chain in model:
        for res in chain:
            if res.id[1] in [202, 315, 363]:
                for atom in res:
                    if atom.name == "CA":
                        active_coords.append(atom.coord)
                        print(f"ACSS2 Active Site Residue {res.resname}{res.id[1]} CA: {atom.coord}")

acss2_center = np.mean(active_coords, axis=0)
acss2_box = {
    "center_x": float(acss2_center[0]),
    "center_y": float(acss2_center[1]),
    "center_z": float(acss2_center[2]),
    "size_x": 22.0,
    "size_y": 22.0,
    "size_z": 22.0
}
print("ACSS2 Docking Grid Box:", acss2_box)

acss2_pdbqt = os.path.join(DATA_DIR, "acss2_clean.pdbqt")
subprocess.run([
    "obabel", "-ipdb", os.path.join(DATA_DIR, "AF-Q9NR19-F1.pdb"),
    "-opdbqt", "-O", acss2_pdbqt,
    "-xr", "-xh"
], check=True)
print(f"Generated {acss2_pdbqt}")

# Save grid configs
import json
with open(os.path.join(DATA_DIR, "grid_boxes.json"), "w") as f:
    json.dump({"ffar2": ffar2_box, "acss2": acss2_box}, f, indent=2)
print("Saved grid_boxes.json successfully.")
