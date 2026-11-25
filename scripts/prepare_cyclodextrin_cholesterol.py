#!/usr/bin/env python3
"""
Generate 3D structures and PDBQT for:
1. beta-Cyclodextrin (host molecule)
2. Cholesterol (guest molecule)
and run host-guest inclusion docking.
"""
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import requests

DATA_DIR = "/run/media/julian/ML6/arteriosklerosis/data"
LIGAND_DIR = os.path.join(DATA_DIR, "ligands")
os.makedirs(LIGAND_DIR, exist_ok=True)

# 1. Download experimental high-resolution structure of beta-Cyclodextrin from PDB (e.g., PDB 3CGT, 1BFN or 4RER ligand BCD)
# Or fetch PubChem 3D SDF for beta-cyclodextrin (CID 444899) and cholesterol (CID 5997)
print("Fetching 3D conformers from PubChem...")
cd_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/444899/SDF?record_type=3d"
chol_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/5997/SDF?record_type=3d"

cd_sdf = os.path.join(DATA_DIR, "beta_cyclodextrin.sdf")
chol_sdf = os.path.join(LIGAND_DIR, "cholesterol.sdf")

r_cd = requests.get(cd_url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30)
if r_cd.status_code == 200:
    with open(cd_sdf, "wb") as f:
        f.write(r_cd.content)
    print(f"Downloaded beta-cyclodextrin 3D SDF ({len(r_cd.content)} bytes)")
else:
    print(f"Failed to fetch cyclodextrin: HTTP {r_cd.status_code}")

r_chol = requests.get(chol_url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30)
if r_chol.status_code == 200:
    with open(chol_sdf, "wb") as f:
        f.write(r_chol.content)
    print(f"Downloaded cholesterol 3D SDF ({len(r_chol.content)} bytes)")
else:
    print(f"Failed to fetch cholesterol: HTTP {r_chol.status_code}")

# Convert both to PDBQT with OpenBabel
cd_pdbqt = os.path.join(DATA_DIR, "beta_cyclodextrin.pdbqt")
chol_pdbqt = os.path.join(LIGAND_DIR, "cholesterol.pdbqt")

subprocess.run(["obabel", "-isdf", cd_sdf, "-opdbqt", "-O", cd_pdbqt, "-xr", "-xh"], check=True)
subprocess.run(["obabel", "-isdf", chol_sdf, "-opdbqt", "-O", chol_pdbqt, "--partialcharge", "gasteiger"], check=True)

# Sanitize receptor PDBQT
with open(cd_pdbqt) as f:
    lines = f.readlines()
clean_lines = [l for l in lines if l.startswith(("ATOM", "HETATM", "TER", "END"))]
with open(cd_pdbqt, "w") as f:
    f.writelines(clean_lines)

print("Cyclodextrin host and cholesterol guest prepared successfully for docking!")
