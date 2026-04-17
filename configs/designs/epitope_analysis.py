import os
import glob
import pymol
from pymol import cmd

# Launch headless
pymol.finish_launching(['pymol', '-cq'])

cif_folder = "/Users/nhinguyen/Desktop/LoveLab/BoltzGen/boltzgen/experiments/designs/Polio_VLP"
cif_files = sorted(glob.glob(os.path.join(cif_folder, "*.cif")))

if not cif_files:
    raise ValueError("No CIF files found.")

reference_view = None
reference_name = None

for i, cif in enumerate(cif_files):
    name = os.path.basename(cif).replace(".cif", "")
    print(f"Processing {name}")

    cmd.load(cif, name)

    cmd.select("fab", f"{name} and chain A")
    cmd.select("vlp", f"{name} and chain B")

    if i == 0:
        reference_name = name

        # Orient so VLP faces viewer
        cmd.orient("vlp")
        cmd.zoom("vlp")
        cmd.turn("y", 180)

        reference_view = cmd.get_view()
    else:
        cmd.align(name, reference_name)
        cmd.set_view(reference_view)

    cmd.hide("everything", name)
    cmd.show("sticks", "vlp")
    cmd.show("surface", "fab")

    cmd.color("orange", "vlp")
    cmd.color("gray", "fab")

    cmd.label("vlp and name CA", f'"{name}"')

# Enable grid view
cmd.set("grid_mode", 1)

# Render combined image
cmd.png("comparison_grid.png", width=2400, height=1800, dpi=300, ray=1)

print("Saved: comparison_grid.png")