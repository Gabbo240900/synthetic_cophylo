

import os
import re

# Base directory containing the Dataset files
BASE_DIR = "/Users/gabriele/synthetic_cophylo/generate_alcala/generate_host_freq/generate_host_freq/generated_trees/Dataset"

# Regex to match Dataset_x.tgl
pattern = re.compile(r"^Dataset_(\d+)\.tgl$")

for filename in os.listdir(BASE_DIR):
    match = pattern.match(filename)
    if not match:
        continue

    old_number = int(match.group(1))
    new_number = old_number + 900

    old_path = os.path.join(BASE_DIR, filename)
    new_filename = f"Dataset_{new_number}.tgl"
    new_path = os.path.join(BASE_DIR, new_filename)

    if os.path.exists(new_path):
        raise FileExistsError(f"Target file already exists: {new_filename}")

    os.rename(old_path, new_path)
    print(f"Renamed {filename} -> {new_filename}")