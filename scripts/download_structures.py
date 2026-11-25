#!/usr/bin/env python3
"""
Download relevant macromolecular structures for molecular modeling of acetate in atherosclerosis.
"""
import os
import requests

DATA_DIR = "/run/media/julian/ML6/arteriosklerosis/data"
os.makedirs(DATA_DIR, exist_ok=True)

URLS = {
    "8J24.pdb": "https://files.rcsb.org/download/8J24.pdb",  # Human FFAR2-Gi with acetate (Cryo-EM 2.6A)
    "AF-Q9NR19-F1.pdb": "https://alphafold.ebi.ac.uk/files/AF-Q9NR19-F1-model_v4.pdb",  # Human ACSS2 AlphaFold
    "4CFE.pdb": "https://files.rcsb.org/download/4CFE.pdb",  # Human AMPK core bound with AMP
    "1HW9.pdb": "https://files.rcsb.org/download/1HW9.pdb",  # Human HMG-CoA Reductase catalytic domain
}

for filename, url in URLS.items():
    filepath = os.path.join(DATA_DIR, filename)
    if not os.path.exists(filepath):
        print(f"Downloading {filename} from {url}...")
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            with open(filepath, "wb") as f:
                f.write(resp.content)
            print(f"Saved {filename} ({len(resp.content)} bytes)")
        else:
            print(f"Failed to download {filename}: HTTP {resp.status_code}")
    else:
        print(f"File {filename} already exists.")
