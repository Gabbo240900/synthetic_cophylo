import os
import re

def parse_tgl_file(filepath):
    host_newick = None
    parasite_newick = None
    associations = []
    sim_stats = {}

    in_host_block = False
    in_parasite_block = False
    in_distribution_block = False

    with open(filepath, 'r') as file:
        lines = file.readlines()

    for line in lines:
        stripped = line.strip()

        # Identify block starts
        if stripped.startswith("BEGIN HOST") or stripped.startswith("# HOST_TREE"):
            in_host_block = True
            continue
        elif stripped.startswith("BEGIN PARASITE") or stripped.startswith("# PARASITE_TREE"):
            in_parasite_block = True
            continue
        elif stripped.startswith("BEGIN DISTRIBUTION") or stripped.startswith("# ASSOCIATIONS"):
            in_distribution_block = True
            continue

        # Identify block ends
        if stripped.startswith("ENDBLOCK"):
            in_host_block = in_parasite_block = in_distribution_block = False
            continue

        # Extract Newick strings
        if in_host_block and stripped and not stripped.startswith("#"):
            host_newick = stripped.strip(" ;")

        elif in_parasite_block and stripped and not stripped.startswith("#"):
            parasite_newick = stripped.strip(" ;")

        # Extract associations
        elif in_distribution_block and not stripped.startswith("#") and " " in stripped:
            parts = stripped.split()
            if len(parts) == 2:
                parasite, host = parts
                associations.append((parasite, host))

    # Parse PARAMETERS block
    param_stats = {}
    in_param_block = False

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("# PARAMETERS"):
            in_param_block = True
            continue
        if stripped.startswith("ENDBLOCK") and in_param_block:
            in_param_block = False
            continue

        if in_param_block and "=" in stripped:
            key, value = map(str.strip, stripped.split("=", 1))
            key = key.replace(" ", "_")
            try:
                val = float(value)
                param_stats[key] = val
            except ValueError:
                param_stats[key] = value

    sim_stats = param_stats

    data = {
        "filename": os.path.basename(filepath),
        "host_newick": host_newick,
        "parasite_newick": parasite_newick,
        "associations": associations,
        "sim_stats": sim_stats
    }

    return data

def parse_all_tgl_in_folder(folder_path):
    datasets = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".tgl"):
            full_path = os.path.join(folder_path, filename)
            try:
                data = parse_tgl_file(full_path)
                datasets.append(data)
            except Exception as e:
                print(f"Error parsing {filename}: {e}")
    return datasets

datasets_coala_cosp = parse_all_tgl_in_folder('/Users/gabriele/synthetic_cophylo/generate_alcala/generated_trees_highCosp/Dataset')
print(datasets_coala_cosp)