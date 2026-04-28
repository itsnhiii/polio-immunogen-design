#!/usr/bin/env python3

import os
import glob
import subprocess

# =========================
# CONFIG
# =========================

INPUT_DIR = "."
OUTPUT_DIR = "."

# 🔥 IMPORTANT: point to working PyMOL binary
PYMOL_BIN = "/Applications/PyMOL.app/Contents/MacOS/PyMOL"

# =========================
# FUNCTION
# =========================

def convert(cif_file):
    name = os.path.basename(cif_file).replace(".cif", "")
    pdb_file = os.path.join(OUTPUT_DIR, f"{name}.pdb")

    if os.path.exists(pdb_file):
        print(f"Skipping (exists): {pdb_file}")
        return

    cmd = f'{PYMOL_BIN} -cq -d "load {cif_file}; save {pdb_file}; quit"'
    subprocess.run(cmd, shell=True, check=True)

    print(f"Converted: {cif_file} → {pdb_file}")


# =========================
# MAIN
# =========================

def main():
    files = glob.glob(os.path.join(INPUT_DIR, "*.cif"))

    print(f"Found {len(files)} CIF files")

    for f in files:
        try:
            convert(f)
        except Exception as e:
            print(f"Failed on {f}: {e}")

    print("Done.")


if __name__ == "__main__":
    main()