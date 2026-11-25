#!/usr/bin/env python3
"""
Prepare ligands (Acetate, Formate, Propionate, Butyrate) with 3D coordinates and PDBQT formatting using RDKit and OpenBabel.
"""
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

OUTPUT_DIR = "/run/media/julian/ML6/arteriosklerosis/data/ligands"
os.makedirs(OUTPUT_DIR, exist_ok=True)

LIGANDS = {
    "acetate": ("CC(=O)[O-]", "Acetate (ionized, physiological)"),
    "acetic_acid_neutral": ("CC(=O)O", "Acetic Acid (neutral)"),
    "formate": ("C(=O)[O-]", "Formate (C1)"),
    "propionate": ("CCC(=O)[O-]", "Propionate (C3)"),
    "butyrate": ("CCCC(=O)[O-]", "Butyrate (C4)"),
}

for name, (smiles, desc) in LIGANDS.items():
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    
    # Save SDF
    sdf_path = os.path.join(OUTPUT_DIR, f"{name}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    
    # Convert to PDBQT with Open Babel
    pdbqt_path = os.path.join(OUTPUT_DIR, f"{name}.pdbqt")
    cmd = [
        "obabel", "-isdf", sdf_path,
        "-opdbqt", "-O", pdbqt_path,
        "--partialcharge", "gasteiger"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0:
        print(f"Generated {pdbqt_path} ({desc})")
    else:
        print(f"Error converting {name}: {res.stderr}")
