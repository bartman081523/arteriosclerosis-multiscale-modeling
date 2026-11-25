#!/usr/bin/env python3
"""
Prepare 3D structures and PDBQT for supermarket bioactive compounds:
1. S-Allylcysteine (Aged Garlic / Schwarzer Knoblauch)
2. Urolithin A (Granatapfel gut metabolite)
3. Hydroxytyrosol (Natives Olivenöl Extra / EVOO)
4. Menaquinone-7 fragment (Vitamin K2 from Natto)
"""
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

DATA_DIR = "/run/media/julian/ML6/arteriosklerosis/data"
LIGAND_DIR = os.path.join(DATA_DIR, "ligands")
os.makedirs(LIGAND_DIR, exist_ok=True)

BIOACTIVES = {
    "s_allylcysteine": ("C=CCSC[C@@H](C(=O)O)N", "S-Allylcysteine (Aged Garlic Extract)"),
    "urolithin_a": ("Oc1ccc2c(c1)oc(=O)c1cc(O)ccc12", "Urolithin A (Pomegranate/Ellagitannin Metabolite)"),
    "hydroxytyrosol": ("c1cc(c(cc1CCO)O)O", "Hydroxytyrosol (Extra Virgin Olive Oil)"),
    "menaquinone_k2": ("CC1=C(C(=O)c2ccccc2C1=O)CC=C(C)C", "Menaquinone core / Vitamin K2 (Natto)"),
    "allicin": ("C=CCSS(=O)CC=C", "Allicin (Fresh Garlic)"),
}

for name, (smiles, desc) in BIOACTIVES.items():
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    
    sdf_path = os.path.join(LIGAND_DIR, f"{name}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    
    pdbqt_path = os.path.join(LIGAND_DIR, f"{name}.pdbqt")
    cmd = ["obabel", "-isdf", sdf_path, "-opdbqt", "-O", pdbqt_path, "--partialcharge", "gasteiger"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0:
        print(f"Generated {pdbqt_path} ({desc})")
    else:
        print(f"Error for {name}: {res.stderr}")
