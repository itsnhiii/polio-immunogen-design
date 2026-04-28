import os
import glob
import re
import pandas as pd
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.align import alignto
from MDAnalysis.analysis.rms import rmsd

# =========================
# CONFIG
# =========================

FOLDER = "./designs_pdb"
NATIVE_FILE = "./native.pdb"

NATIVE_KEY_SEL = "segid E or segid F"
NATIVE_LOCK_SEL = "segid A or segid B or segid C"

DESIGN_KEY_SEL = "segid A or segid B"
DESIGN_LOCK_SEL = "segid C"

CONTACT_CUTOFF = 4.0

CLUSTER_A = [101, 102, 103, 104, 105, 107, 108, 109, 110]
CLUSTER_B = [25, 26, 30, 54]
CLUSTER_C = [94, 95]

CSV_OUT = os.path.join(FOLDER, "analysis_results.csv")

# =========================
# HELPERS
# =========================

def natural_key(s):
    """Natural sort key: run2 before run10."""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]

def get_contacts(key_atoms, lock_atoms):
    contacting = set()
    for res in key_atoms.residues:
        for atom_i in res.atoms:
            dists = np.linalg.norm(lock_atoms.positions - atom_i.position, axis=1)
            if np.any(dists < CONTACT_CUTOFF):
                contacting.add(res.resid)
                break
    return contacting

def jaccard(a, b):
    return len(a & b) / len(a | b) if len(a | b) > 0 else 0.0

def cluster_coverage(contact_set, cluster):
    return len(set(cluster) & contact_set) / len(cluster)

def min_dist(key_atoms, lock_atoms, cluster):
    atoms = key_atoms.select_atoms(" or ".join([f"resid {r}" for r in cluster]))
    if len(atoms) == 0 or len(lock_atoms) == 0:
        return np.nan
    dmin = np.inf
    for atom in atoms:
        dists = np.linalg.norm(lock_atoms.positions - atom.position, axis=1)
        dmin = min(dmin, np.min(dists))
    return dmin

def missing_hotspots(contact_set, cluster):
    return sorted(set(cluster) - set(contact_set))

# =========================
# ANALYSIS
# =========================

def analyze(native_u, pdb_file):
    u_des = mda.Universe(pdb_file)

    key_nat = native_u.select_atoms(NATIVE_KEY_SEL)
    lock_nat = native_u.select_atoms(NATIVE_LOCK_SEL)

    key_des = u_des.select_atoms(DESIGN_KEY_SEL)
    lock_des = u_des.select_atoms(DESIGN_LOCK_SEL)

    # Align design onto native using antibody only.
    # match_atoms=False avoids exact atom-count matching requirements.
    alignto(
        u_des,
        native_u,
        select=(DESIGN_KEY_SEL, NATIVE_KEY_SEL),
        match_atoms=False
    )

    nat_contacts = get_contacts(key_nat, lock_nat)
    des_contacts = get_contacts(key_des, lock_des)

    jacc = jaccard(nat_contacts, des_contacts)

    # Global lock RMSD after key alignment.
    # This is a rough metric only since native lock and design lock are not atom-matched objects.
    n = min(len(lock_nat.positions), len(lock_des.positions))
    lock_rmsd = rmsd(lock_des.positions[:n], lock_nat.positions[:n])

    A_cov = cluster_coverage(des_contacts, CLUSTER_A)
    B_cov = cluster_coverage(des_contacts, CLUSTER_B)
    C_cov = cluster_coverage(des_contacts, CLUSTER_C)

    A_d = min_dist(key_des, lock_des, CLUSTER_A)
    B_d = min_dist(key_des, lock_des, CLUSTER_B)
    C_d = min_dist(key_des, lock_des, CLUSTER_C)

    A_missing = missing_hotspots(des_contacts, CLUSTER_A)
    B_missing = missing_hotspots(des_contacts, CLUSTER_B)
    C_missing = missing_hotspots(des_contacts, CLUSTER_C)

    return {
        "file": os.path.basename(pdb_file),
        "jaccard": jacc,
        "rmsd": lock_rmsd,
        "A_cov": A_cov,
        "B_cov": B_cov,
        "C_cov": C_cov,
        "A_dist": A_d,
        "B_dist": B_d,
        "C_dist": C_d,
        "A_missing": ",".join(map(str, A_missing)),
        "B_missing": ",".join(map(str, B_missing)),
        "C_missing": ",".join(map(str, C_missing)),
    }

# =========================
# PLOTTING
# =========================

def save_coverage_plot(df, outpath):
    dfp = df.sort_values("file", key=lambda col: col.map(natural_key)).reset_index(drop=True)
    x = np.arange(len(dfp))
    labels = dfp["file"]

    fig, axes = plt.subplots(3, 1, figsize=(max(10, len(dfp) * 0.35), 9), sharex=True, sharey=True)

    axes[0].plot(x, dfp["A_cov"], marker="o")
    axes[0].set_ylabel("A coverage")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Coverage by run (increasing run order)")

    axes[1].plot(x, dfp["B_cov"], marker="o")
    axes[1].set_ylabel("B coverage")
    axes[1].set_ylim(0, 1.05)

    axes[2].plot(x, dfp["C_cov"], marker="o")
    axes[2].set_ylabel("C coverage")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=90)
    axes[2].set_xlabel("Run")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

def save_distance_plot(df, outpath):
    dfp = df.sort_values("file", key=lambda col: col.map(natural_key), ascending=False).reset_index(drop=True)
    x = np.arange(len(dfp))
    labels = dfp["file"]

    max_y = np.nanmax([dfp["A_dist"].max(), dfp["B_dist"].max(), dfp["C_dist"].max()])

    fig, axes = plt.subplots(3, 1, figsize=(max(10, len(dfp) * 0.35), 9), sharex=True, sharey=True)

    axes[0].plot(x, dfp["A_dist"], marker="o")
    axes[0].set_ylabel("A min dist (Å)")
    axes[0].set_ylim(0, max_y * 1.05)
    axes[0].set_title("Minimum distance by run (decreasing run order)")

    axes[1].plot(x, dfp["B_dist"], marker="o")
    axes[1].set_ylabel("B min dist (Å)")
    axes[1].set_ylim(0, max_y * 1.05)

    axes[2].plot(x, dfp["C_dist"], marker="o")
    axes[2].set_ylabel("C min dist (Å)")
    axes[2].set_ylim(0, max_y * 1.05)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=90)
    axes[2].set_xlabel("Run")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

def save_single_metric_plot(df, metric, ylabel, title, outpath):
    dfp = df.sort_values("file", key=lambda col: col.map(natural_key)).reset_index(drop=True)
    x = np.arange(len(dfp))
    labels = dfp["file"]

    plt.figure(figsize=(max(10, len(dfp) * 0.35), 4.5))
    plt.plot(x, dfp[metric], marker="o")
    plt.xticks(x, labels, rotation=90)
    plt.xlabel("Run")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

def save_score_hist(df, outpath):
    plt.figure(figsize=(6, 4.5))
    plt.hist(df["score"], bins=20)
    plt.xlabel("Score")
    plt.ylabel("Count")
    plt.title("Design score distribution")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()

# =========================
# MAIN
# =========================

def main():
    native_u = mda.Universe(NATIVE_FILE)

    pdbs = sorted(glob.glob(os.path.join(FOLDER, "*.pdb")), key=natural_key)
    results = []

    for f in pdbs:
        base = os.path.basename(f).lower()
        if base == os.path.basename(NATIVE_FILE).lower():
            continue
        print(f"Processing {f}")
        res = analyze(native_u, f)
        results.append(res)

    df = pd.DataFrame(results)

    # Composite score
    df["score"] = (
        2.0 * df["A_cov"]
        + 1.5 * df["B_cov"]
        + 1.5 * df["C_cov"]
        - 0.1 * (df["B_dist"] + df["C_dist"])
    )

    # Save CSV
    df.to_csv(CSV_OUT, index=False)
    print(f"Saved {CSV_OUT}")

    # Save plots
    save_coverage_plot(df, os.path.join(FOLDER, "coverage_subplots.png"))
    save_distance_plot(df, os.path.join(FOLDER, "distance_subplots.png"))
    save_single_metric_plot(
        df,
        metric="jaccard",
        ylabel="Jaccard",
        title="Jaccard by run",
        outpath=os.path.join(FOLDER, "jaccard_by_run.png"),
    )
    save_single_metric_plot(
        df,
        metric="rmsd",
        ylabel="RMSD (Å)",
        title="RMSD by run",
        outpath=os.path.join(FOLDER, "rmsd_by_run.png"),
    )
    save_score_hist(df, os.path.join(FOLDER, "score_hist.png"))

    print("Plots saved.")

if __name__ == "__main__":
    main()