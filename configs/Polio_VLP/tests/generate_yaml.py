import os

# =========================
# CONFIG
# =========================

# Constraint sets
binding_sets = {
    # "A": {
    #     "E": "100..103,105..109"
    # },
    # "B": {
    #     "E": "97..103,105..109"
    # },
    # "C": {
    #     "E": "97..103,105..109",
    #     "F": "93..95"
    # },
    # "D": {
    #     "E": "29,52,97..103,105..109",
    #     "F": "93..95"
    # },
    # "E": {
    #     "E": "100,105,108"
    # },
    "F": {
        "E": "101,103,106",
        "F": "54,95"    
    }
}

# Length ranges
lengths = [
    (80, 110),
    (100, 140),
    (120, 160),
    (160, 200)
]

# =========================
# FUNCTION
# =========================

def make_yaml(binding_config, length_range):
    binding_block = ""

    for chain_id, residues in binding_config.items():
        binding_block += f"""
        - chain:
            id: {chain_id}
            binding: {residues}
"""

    yaml_str = f"""
entities:
  - file:
      path: 8E8Z.cif

      include:
        - chain:
            id: E
        - chain:
            id: F

      binding_types:{binding_block}

      structure_groups:
        - group:
            visibility: 1
            id: E
        - group:
            visibility: 1
            id: F
            
  - protein:
      id: P
      sequence: {length_range[0]}..{length_range[1]}
"""
    return yaml_str.strip()

# =========================
# GENERATE FILES
# =========================

for set_name, binding in binding_sets.items():
    for lmin, lmax in lengths:

        #filename = f"8e8z_{set_name}_{lmin}_{lmax}.yaml" # ending yaml file name with number is illegal
        filename = f"8e8z_{set_name}_len{lmin}to{lmax}.yaml"
        yaml_content = make_yaml(binding, (lmin, lmax))

        with open(filename, "w") as f:
            f.write(yaml_content)

        print(f"Created {filename}")

print("\nDone.")