#!/usr/bin/env python3
"""
Convert PyMOL auth-style residues (chain auth_resi resn)
to BoltzGen label_seq_id YAML using mmCIF.

Example input line:
H 81 GLN
L 45 ILE
"""

import gemmi
from collections import defaultdict

### ========= USER INPUT ========= ###

CIF_FILE = "8E8Z.cif"   # <-- change if needed

# Paste PyMOL output here (auth_asym_id, auth_seq_id, resname)
PYMOL_RESIDUES = """
H 81 GLN
H 119 SER
H 120 ARG
H 121 PRO
H 122 SER
H 123 ALA
H 124 TRP
H 125 VAL
H 127 ARG
H 128 SER
H 129 LEU
H 130 TYR
L 45 ILE
L 49 GLY
L 50 TYR
L 73 THR
L 74 ASN
L 75 ARG
L 90 ASN
L 113 ARG
L 114 SER
L 115 ASN
"""

### =============================== ###


def parse_pymol_input(text):
    residues = []
    for line in text.strip().splitlines():
        chain, resi, resn = line.split()
        residues.append((chain, int(resi)))
    return residues


def build_auth_to_label_map(block):
    """
    Build mapping:
    (auth_asym_id, auth_seq_id) -> (label_asym_id, label_seq_id)
    """
    mapping = {}

    loop = block.find_loop("_pdbx_poly_seq_scheme.auth_seq_id")
    tags = loop.tags

    idx = {tag: i for i, tag in enumerate(tags)}

    for row in loop:
        auth_chain = row[idx["_pdbx_poly_seq_scheme.auth_asym_id"]]
        auth_seq   = int(row[idx["_pdbx_poly_seq_scheme.auth_seq_id"]])
        label_chain = row[idx["_pdbx_poly_seq_scheme.asym_id"]]
        label_seq   = int(row[idx["_pdbx_poly_seq_scheme.seq_id"]])

        mapping[(auth_chain, auth_seq)] = (label_chain, label_seq)

    return mapping


def main():
    pymol_residues = parse_pymol_input(PYMOL_RESIDUES)

    doc = gemmi.cif.read_file(CIF_FILE)
    block = doc.sole_block()

    auth_to_label = build_auth_to_label_map(block)

    output = defaultdict(set)
    missing = []

    for auth_chain, auth_resi in pymol_residues:
        key = (auth_chain, auth_resi)
        if key not in auth_to_label:
            missing.append(key)
            continue

        label_chain, label_resi = auth_to_label[key]
        output[label_chain].add(label_resi)

    print("\n### BoltzGen YAML (label_seq_id) ###\n")

    for chain in sorted(output):
        residues = sorted(output[chain])
        print(f"- chain:")
        print(f"    id: {chain}")
        print(f"    binding: {','.join(map(str, residues))}\n")

    if missing:
        print("### WARNING: Unmapped residues (likely non-polymer or wrong structure) ###")
        for c, r in missing:
            print(f"  auth_chain={c} auth_resi={r}")


if __name__ == "__main__":
    main()