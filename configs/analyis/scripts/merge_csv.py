import pandas as pd
import glob
import os

exp_name = "Polio_VLP"
analysis_folder = "/Users/nhinguyen/Desktop/LoveLab/BoltzGen/boltzgen/experiments/analyis/"
folder = analysis_folder + exp_name
files = sorted(glob.glob(os.path.join(folder, "*.csv")))

dfs = []

for f in files:
    df = pd.read_csv(f)
    df.insert(1, "source_file", os.path.basename(f))
    dfs.append(df)

merged = pd.concat(dfs, ignore_index=True)

output_file = os.path.join(analysis_folder, exp_name + ".csv")
merged.to_csv(output_file, index=False)

print(f"Merged {len(files)} files into {output_file}")