import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.align import alignto
from MDAnalysis.analysis.rms import rmsd

CONTACT_CUTOFF = 4.0  # Å

# native structure
NATIVE_FILE = "native.pdb"
NATIVE_KEY_SEL = "segid E or segid F"
NATIVE_LOCK_SEL = "segid A or segid B or segid C"

# design structure
DESIGN_FILE = "A_120-160.pdb"
DESIGN_KEY_SEL = "segid A or segid B"
DESIGN_LOCK_SEL = "segid C"

# key hotspot clusters on the antibody side
# native numbering after shift convention (-22 HC, -20 HC)
# CLUSTER_A = [101, 103, 106]   # heavy-chain hotspot
CLUSTER_A = [101, 102, 103, 104, 105, 107, 108, 109, 110]     # whole patch
CLUSTER_B = [25, 26, 30, 54]          # light-chain patch 1
CLUSTER_C = [94, 95]              # light-chain patch 2


def get_contacting_key_residues(key_atoms, lock_atoms, cutoff=4.5):
    contacting = set()
    for res in key_atoms.residues:
        for atom_i in res.atoms:
            dists = np.linalg.norm(lock_atoms.positions - atom_i.position, axis=1)
            if np.any(dists < cutoff):
                contacting.add(res.resid)
                break
    return contacting


def jaccard(a, b):
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def cluster_coverage(contact_set, cluster):
    cluster = set(cluster)
    if not cluster:
        return 0.0
    return len(contact_set & cluster) / len(cluster)


def min_distance_to_cluster(key_atoms, lock_atoms, cluster_resids):
    cluster_atoms = key_atoms.select_atoms(" or ".join([f"resid {r}" for r in cluster_resids]))
    if len(cluster_atoms) == 0 or len(lock_atoms) == 0:
        return np.nan
    dmin = np.inf
    for atom in cluster_atoms:
        dists = np.linalg.norm(lock_atoms.positions - atom.position, axis=1)
        dmin = min(dmin, np.min(dists))
    return dmin


def main():
    u_nat = mda.Universe(NATIVE_FILE)
    u_des = mda.Universe(DESIGN_FILE)

    key_nat = u_nat.select_atoms(NATIVE_KEY_SEL)
    key_des = u_des.select_atoms(DESIGN_KEY_SEL)

    lock_nat = u_nat.select_atoms(NATIVE_LOCK_SEL)
    lock_des = u_des.select_atoms(DESIGN_LOCK_SEL)

    # align design onto native using antibody/key only
    alignto(u_des, u_nat, select=(DESIGN_KEY_SEL, NATIVE_KEY_SEL))

    # key-side contact sets
    nat_contacts = get_contacting_key_residues(key_nat, lock_nat, CONTACT_CUTOFF)
    des_contacts = get_contacting_key_residues(key_des, lock_des, CONTACT_CUTOFF)

    # global lock RMSD after key alignment
    # note: native lock is the real antigen, design lock is your designed immunogen,
    # so this is only a rough shape metric, not a strict residue-to-residue equivalence metric
    n = min(len(lock_nat.positions), len(lock_des.positions))
    lock_rmsd = rmsd(lock_des.positions[:n], lock_nat.positions[:n])

    # key-centric similarity
    jacc = jaccard(nat_contacts, des_contacts)

    # hotspot coverage on design
    A_cov = cluster_coverage(des_contacts, CLUSTER_A)
    B_cov = cluster_coverage(des_contacts, CLUSTER_B)
    C_cov = cluster_coverage(des_contacts, CLUSTER_C)

    # minimum distances from design lock to each hotspot cluster
    A_dmin = min_distance_to_cluster(key_des, lock_des, CLUSTER_A)
    B_dmin = min_distance_to_cluster(key_des, lock_des, CLUSTER_B)
    C_dmin = min_distance_to_cluster(key_des, lock_des, CLUSTER_C)

    print("=== Design vs Native Recognition Metrics ===")
    print(f"Jaccard contact similarity (key-centric): {jacc:.3f}")
    print(f"Lock RMSD after key alignment: {lock_rmsd:.2f} Å")
    print("")
    print("Cluster coverage:")
    print(f"  A coverage: {A_cov:.2f}")
    print(f"  B coverage: {B_cov:.2f}")
    print(f"  C coverage: {C_cov:.2f}")
    print("")
    print("Minimum distance from design lock to cluster:")
    print(f"  A min distance: {A_dmin:.2f} Å")
    print(f"  B min distance: {B_dmin:.2f} Å")
    print(f"  C min distance: {C_dmin:.2f} Å")


if __name__ == "__main__":
    main()