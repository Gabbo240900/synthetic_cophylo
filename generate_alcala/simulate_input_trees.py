#ALCALA GENERATOR - add time for host conditioning birth death 
import os
import random
import argparse
import pandas as pd
import re
import subprocess
import shutil
from ete3 import Tree
from dendropy.simulate import treesim
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np
from io import StringIO
from Bio import Phylo
from collections import defaultdict
# Track datasets with errors
error_datasets = []


# Clean existing folders if present
for folder in ["alcala_trees", "Dataset", "host_trees", "parasite_trees"]:
    full_path = os.path.join("generated_trees", folder)
    if os.path.exists(full_path):
        shutil.rmtree(full_path)

class GenerateHostTree:
    """Class to generate host trees and their corresponding event frequency files."""
    
    def __init__(self, num_trees, min_leaves, max_leaves, output_dir, tmrca=2):
        self.num_trees = num_trees  
        self.min_leaves = min_leaves  
        self.max_leaves = max_leaves  
        self.output_dir = output_dir  
        self.tmrca = tmrca

        # Define host tree and frequency directories
        self.host_tree_dir = os.path.join("generated_trees", "host_trees")
        self.alcala_tree_dir = os.path.join("generated_trees", "alcala_trees")
        os.makedirs(self.alcala_tree_dir, exist_ok=True)

        # Ensure output directories exist
        os.makedirs(self.host_tree_dir, exist_ok=True)
        
    def run_cophylo_simulation(self, host_tree_index):
        base_path = os.path.join(self.alcala_tree_dir, f"host_tree_{host_tree_index}")
        input_prefix = os.path.join(base_path, "")
        # Define parasite output folder and prefix
        parasite_folder = os.path.join(self.alcala_tree_dir, f"parasite_tree_{host_tree_index}")
        os.makedirs(parasite_folder, exist_ok=True)
        output_prefix = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}")
        distrib_path = os.path.join(base_path, f"host_tree_{host_tree_index}_distrib.txt")
        s = random.uniform(0.2, 0.4)
        c = random.uniform(0.2, 0.4)

        cmd = [
            "software/bin/cophylo.out",
            "-l", "0.7",
            "-m", "0.6",
            '-c', str(c),
            '-s', str(s),
            "-t", str(self.tmrca),
            "-i", input_prefix,
            "-o", output_prefix,
            "-N", "1",
            "-S", str(random.randint(1, 1000000)),
            "-P", distrib_path
        ]
        try:
            subprocess.run(cmd, check=True)
            edge_txt_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}edges_0.txt")
            length_txt_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}edgelength_0.txt")
            parasite_newick_dir = os.path.join("generated_trees", "parasite_trees")
            os.makedirs(parasite_newick_dir, exist_ok=True)
            newick_output_path = os.path.join(parasite_newick_dir, f"parasite_tree_{host_tree_index}.nwk")
            self.txt_to_newick(edge_txt_path, length_txt_path, newick_output_path)
            host_nwk_path = os.path.join(self.host_tree_dir, f"host_tree_{host_tree_index}.nwk")
            hpassoc_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}hpassoc_0.txt")
            params_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}params_0.txt")
            dataset_dir = os.path.join("generated_trees", "Dataset")
            os.makedirs(dataset_dir, exist_ok=True)
            dataset_path = os.path.join(dataset_dir, f"Dataset_{host_tree_index}.tgl")
            self.create_tgl_dataset(host_nwk_path, newick_output_path, hpassoc_path, params_path, dataset_path)
        except subprocess.CalledProcessError as e:
            print("Simulation failed for tree", host_tree_index)
            print("Command:", ' '.join(cmd))
            print("Error:", e)


    def generate_random_tree(self, num_leaves, prefix="H"):
        s = te.species_tree_n_age(n=num_leaves, model='BDP', age=self.tmrca, birth_rate=0.7, death_rate=0.45)
        s_nwk = to_newick(s). strip()
        if not s_nwk.endswith(";"):
            s_nwk += ";"
        s_nwk = re.sub(r'(?<=[(,)])(?!(H))(\d+):', r'H\2:', s_nwk)
        tree = Tree(s_nwk, format=1)
        for idx, leaf in enumerate(tree.iter_leaves()):
            leaf.name = f"{prefix}{idx+1}"

        # Assign valid names to internal nodes if missing or empty
        for node in tree.traverse("postorder"):
            if not node.is_leaf() and (not node.name or not node.name.strip()):
                node.name = f"IN{random.randint(10000, 99999)}"
        return tree

    def save_tree(self, tree, filename):
        dir_branches = self.host_tree_dir
        os.makedirs(dir_branches, exist_ok=True)
        tree_path2 = os.path.join(dir_branches, filename)
        tree.write(outfile=tree_path2, format=5)
        
    def newick_to_txt(self, newick_path, edge_txt_path, length_txt_path):
        tree = Tree(newick_path, format=5)
        edges = []
        lengths = []

        # Assign consistent integer names to internal nodes
        counter = 1
        num_leaves = len(tree.get_leaves())
        for node in tree.traverse("postorder"):
            if not node.name or not node.name.strip():
                node.name = str(num_leaves + counter)
                counter += 1

        for node in tree.traverse("postorder"):
            if not node.is_root():
                parent = node.up
                # Remove 'H' prefix from both parent and node names
                edges.append(f"{parent.name.replace('H', '')} {node.name.replace('H', '')}")
                # Ensure all lengths are positive and not None
                if node.dist is None or node.dist <= 0.0:
                    node.dist = 0.01
                # Check for invalid branch lengths
                if np.isnan(node.dist) or node.dist > 1e6:
                    print(f"WARNING: Skipping invalid length {node.dist} for edge {parent.name} -> {node.name}")
                    continue
                lengths.append(str(node.dist))

        with open(edge_txt_path, "w") as edge_file:
            edge_file.write("\n".join(edges))

        with open(length_txt_path, "w") as length_file:
            length_file.write("\n".join(lengths))

    def txt_to_newick(self, edge_txt_path, length_txt_path, output_newick_path):

        with open(edge_txt_path) as f:
            edges = [line.strip().split() for line in f.readlines()]
        with open(length_txt_path) as f:
            lengths = [float(x.strip()) for x in f.readlines()]

        child_map = defaultdict(list)
        parent_map = {}
        nodes = {}

        for i, (parent, child) in enumerate(edges):
            child_map[parent].append((child, lengths[i]))
            parent_map[child] = parent
            nodes[parent] = nodes.get(parent, Tree(name=parent))
            nodes[child] = Tree(name=child)
            nodes[child].dist = lengths[i]

        all_nodes = set(nodes.keys())
        roots = [n for n in all_nodes if n not in parent_map]
        root_name = roots[0]
        root = nodes[root_name]

        # Prefix all node names with "P"
        for node in nodes.values():
            node.name = f"P{node.name}"

        def attach_children(node_name):
            for child_name, _ in child_map.get(node_name, []):
                nodes[node_name].add_child(nodes[child_name])
                attach_children(child_name)

        attach_children(root_name)
        os.makedirs(os.path.dirname(output_newick_path), exist_ok=True)
        root.write(outfile=output_newick_path, format=5)

    def create_tgl_dataset(self, host_newick_path, parasite_newick_path, hpassoc_path, params_path, output_path):
        with open(host_newick_path) as f:
            host_nwk = f.read().strip()
        with open(parasite_newick_path) as f:
            parasite_nwk = f.read().strip()
        with open(hpassoc_path) as f:
            hpassoc = f.read().strip()
        with open(params_path) as f:
            lines = f.readlines()
            cospeciation = lines[0].strip().split()[0]
            host_switch = lines[0].strip().split()[1]

        with open(output_path, "w") as f:
            f.write("# HOST_TREE\n")
            f.write(host_nwk + "\n")
            f.write("ENDBLOCK\n\n")
            f.write("# PARASITE_TREE\n")
            f.write(parasite_nwk + "\n")
            f.write("ENDBLOCK\n\n")
            f.write("# ASSOCIATIONS\n")

            # Reindex parasite IDs to start from 1
            host_tree = Tree(host_newick_path)
            parasite_ids = []
            for line in hpassoc.splitlines():
                if not line.strip():
                    continue
                pid, _ = map(int, line.strip().split())
                parasite_ids.append(pid)
            unique_parasite_ids = sorted(set(parasite_ids))
            parasite_id_map = {pid: idx + 1 for idx, pid in enumerate(unique_parasite_ids)}

            for line in hpassoc.splitlines():
                if not line.strip():
                    continue
                pid, hid = map(int, line.strip().split())
                pname = f"P{parasite_id_map[pid]}"
                hname = f"H{hid}"
                if hname not in [leaf.name for leaf in host_tree.iter_leaves()]:
                    if output_path not in error_datasets:
                        error_datasets.append(output_path)
                f.write(f"{pname} {hname}\n")

            f.write("ENDBLOCK\n\n")
            f.write("# PARAMETERS\n")
            f.write(f"Cospeciation = {cospeciation}\n")
            f.write(f"Host_switch = {host_switch}\n")
            f.write("ENDBLOCK\n")

    # zipf distribution for specialist parasites
    def save_specialist_distribution(self, tree, folder_path, tree_index): 
        k = len(tree.get_leaves())
        s = 1.6
        ranks = np.arange(1, k + 1, dtype=np.float64)
        probs = 1.0 / np.power(ranks, s)
        probs = probs / probs.sum()

        values = [f"{p:.6f}" for p in probs]
        output_path = os.path.join(folder_path, f"host_tree_{tree_index}_distrib.txt")
        with open(output_path, "w") as f:
            f.write(" ".join(values))

    def generate_and_save_trees(self):


        trees = []
        for i in range(self.num_trees):
            num_leaves = random.randint(self.min_leaves, self.max_leaves)
            tree = self.generate_random_tree(num_leaves)
            tree_filename = f"host_tree_{i + 1}.nwk"
            self.save_tree(tree, tree_filename)
            trees.append((tree, tree_filename, i + 1))
        
        for tree, tree_filename, index in trees:
            tree_folder = os.path.join(self.alcala_tree_dir, f"host_tree_{index}")
            os.makedirs(tree_folder, exist_ok=True)
            newick_path = os.path.join(self.host_tree_dir, tree_filename)
            edge_txt_path = os.path.join(tree_folder, "host_edge.txt")
            length_txt_path = os.path.join(tree_folder, "host_edgelength.txt")
            self.newick_to_txt(newick_path, edge_txt_path, length_txt_path)
            print(f"Converted {tree_filename} to text files in {tree_folder}")

        for tree, _, index in trees:
            tree_folder = os.path.join(self.alcala_tree_dir, f"host_tree_{index}")
            self.save_specialist_distribution(tree, tree_folder, index)
            self.run_cophylo_simulation(index)

        print(f"Generated {self.num_trees} host trees in '{self.host_tree_dir}'.")
        print(f"Converted trees saved in '{self.alcala_tree_dir}'.")


    def get_directories(self):
        """Returns the host tree and frequency directories."""
        return self.host_tree_dir



if __name__ == "__main__":
    import time
    parser = argparse.ArgumentParser(description="Generate host and parasite trees")
    parser.add_argument("--num_trees", type=int, required=True, help="Number of trees to generate")
    parser.add_argument("--min_leaves", type=int, default=15, help="Minimum number of leaves per tree")
    parser.add_argument("--max_leaves", type=int, default=30, help="Maximum number of leaves per tree")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for host trees and frequencies")
    args = parser.parse_args()

    total_start = time.time()

    # Generate Host Trees
    host_generator = GenerateHostTree(args.num_trees, args.min_leaves, args.max_leaves, args.output_dir)
    start = time.time()
    if not os.listdir(os.path.join(args.output_dir, "host_trees")):
        host_generator.generate_and_save_trees()
    else:
        print("Host trees already exist, skipping generation.")
    print(f"Host tree generation took {time.time() - start:.2f} seconds")

    host_dir = host_generator.get_directories()

    print(f"Total execution time: {time.time() - total_start:.2f} seconds")
    # Print errors list
    if error_datasets:
        print("\nDatasets with errors:")
        for err in error_datasets:
            print(f" - {err}")
    else:
        print("\nNo dataset errors encountered.")