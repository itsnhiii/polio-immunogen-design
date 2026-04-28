import MDAnalysis as mda

u = mda.Universe("A_120-160.pdb")

for chain in u.segments:
    atoms = u.select_atoms(f"segid {chain.segid}")
    print(chain.segid, len(atoms.residues))