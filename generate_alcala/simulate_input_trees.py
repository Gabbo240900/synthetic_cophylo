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

# The cophylo binary lives next to this script, so the generator can be run
# from any working directory.
COPHYLO_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "software", "bin", "cophylo.out")


# Birth rates are the same in every regime.
HOST_BIRTH_RATE = 0.7        # lambda_H
SYMBIONT_BIRTH_RATE = 0.7    # lambda_S, cophylo -l

# Per-regime parameters.  "cosp" is the cospeciation probability c (cophylo -c),
# "switch" the host-switch rate s (cophylo -s), both drawn uniformly from their
# range; "host_death" is mu_H (host tree simulation) and "sym_death" is mu_S
# (cophylo -m).
REGIMES = {
    "medium":     {"cosp": (0.2, 0.4),  "switch": (0.2, 0.4),
                   "host_death": 0.45, "sym_death": 0.6},
    "highCosp":   {"cosp": (0.7, 0.9),  "switch": (0.0, 0.05),
                   "host_death": 0.63, "sym_death": 0.645},
    "highSwitch": {"cosp": (0.05, 0.1), "switch": (0.5, 0.7),
                   "host_death": 0.24, "sym_death": 0.66},
}


def clean_output_dirs(base_dir):
    """Remove previously generated folders under base_dir."""
    for folder in ["alcala_trees", "Dataset", "host_trees", "parasite_trees"]:
        full_path = os.path.join(base_dir, folder)
        if os.path.exists(full_path):
            shutil.rmtree(full_path)

class GenerateHostTree:
    """Class to generate host trees and their corresponding event frequency files."""
    
    def __init__(self, num_trees, min_leaves, max_leaves, output_dir, tmrca=2,
                 cosp_range=(0.2, 0.4), switch_range=(0.2, 0.4),
                 host_death=0.45, sym_death=0.6,
                 min_parasites=5, host_distrib="zipf", zipf_s=1.6):
        self.num_trees = num_trees  
        self.min_leaves = min_leaves  
        self.max_leaves = max_leaves  
        self.output_dir = output_dir  
        self.tmrca = tmrca
        self.cosp_range = cosp_range
        self.switch_range = switch_range
        self.host_death = host_death
        self.sym_death = sym_death
        self.min_parasites = min_parasites
        self.host_distrib = host_distrib
        self.zipf_s = zipf_s

        # Define host tree and frequency directories
        self.host_tree_dir = os.path.join(self.output_dir, "host_trees")
        self.alcala_tree_dir = os.path.join(self.output_dir, "alcala_trees")
        self.dataset_dir = os.path.join(self.output_dir, "Dataset")
        self.parasite_tree_dir = os.path.join(self.output_dir, "parasite_trees")
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
        s = random.uniform(*self.switch_range)
        c = random.uniform(*self.cosp_range)

        # Cophylo builds its file names with sprintf into char[100] buffers, so a
        # long -i/-o prefix overflows them and the binary dies without output.
        # Run it from alcala_trees/ and hand it short relative prefixes instead,
        # which keeps the length independent of where the repository lives.
        rel_input_prefix = os.path.join(f"host_tree_{host_tree_index}", "")
        rel_output_prefix = os.path.join(f"parasite_tree_{host_tree_index}",
                                         f"parasite_tree_{host_tree_index}")
        rel_distrib_path = os.path.join(f"host_tree_{host_tree_index}",
                                        f"host_tree_{host_tree_index}_distrib.txt")

        cmd = [
            COPHYLO_BIN,
            "-l", str(SYMBIONT_BIRTH_RATE),
            "-m", str(self.sym_death),
            '-c', str(c),
            '-s', str(s),
            "-t", str(self.tmrca),
            "-i", rel_input_prefix,
            "-o", rel_output_prefix,
            # -N is the minimum number of surviving parasite tips; cophylo
            # retries the simulation until it is reached.
            "-N", str(self.min_parasites),
            "-S", str(random.randint(1, 1000000)),
            "-P", rel_distrib_path
        ]
        try:
            subprocess.run(cmd, check=True, cwd=self.alcala_tree_dir)
            edge_txt_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}edges_0.txt")
            length_txt_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}edgelength_0.txt")
            parasite_newick_dir = self.parasite_tree_dir
            os.makedirs(parasite_newick_dir, exist_ok=True)
            newick_output_path = os.path.join(parasite_newick_dir, f"parasite_tree_{host_tree_index}.nwk")
            self.txt_to_newick(edge_txt_path, length_txt_path, newick_output_path)
            host_nwk_path = os.path.join(self.host_tree_dir, f"host_tree_{host_tree_index}.nwk")
            hpassoc_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}hpassoc_0.txt")
            params_path = os.path.join(parasite_folder, f"parasite_tree_{host_tree_index}params_0.txt")
            dataset_dir = self.dataset_dir
            os.makedirs(dataset_dir, exist_ok=True)
            dataset_path = os.path.join(dataset_dir, f"Dataset_{host_tree_index}.tgl")
            self.create_tgl_dataset(host_nwk_path, newick_output_path, hpassoc_path, params_path, dataset_path)
        except subprocess.CalledProcessError as e:
            print("Simulation failed for tree", host_tree_index)
            print("Command:", ' '.join(cmd))
            print("Error:", e)


    def generate_random_tree(self, num_leaves, prefix="H"):
        s = te.species_tree_n_age(n=num_leaves, model='BDP', age=self.tmrca, birth_rate=HOST_BIRTH_RATE, death_rate=self.host_death)
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
        """Write the host tree as a Cophylo edge/edge-length file pair.

        Cophylo reads the very first integer of host_edge.txt, subtracts one and
        uses that as the number of hosts; it also treats that same value as the
        root node id.  Both only hold under the R `phylo` convention:
        leaves = 1..n, root = n+1, other internal nodes = n+2..2n-1, with the
        edges listed in preorder so that the first edge starts at the root.
        """
        tree = Tree(newick_path, format=5)

        # asymmetree wraps the crown group in an extra branch, which ete3 reads
        # as a root with a single child.  Drop it so the tree is rooted binary.
        while len(tree.children) == 1:
            tree = tree.children[0]
            tree.up = None
            tree.dist = 0.0

        leaves = list(tree.iter_leaves())
        num_leaves = len(leaves)

        for node in tree.traverse():
            if not node.is_leaf() and len(node.children) != 2:
                raise ValueError(
                    f"Host tree {newick_path} is not binary: node {node.name!r} "
                    f"has {len(node.children)} children"
                )

        # Leaves are named H1..Hn upstream; keep those ids so that the host ids
        # Cophylo reports in the association file map back to the leaf names.
        expected_leaf_ids = set(range(1, num_leaves + 1))
        observed_leaf_ids = set()
        for leaf in leaves:
            match = re.fullmatch(r"H(\d+)", leaf.name)
            if not match:
                raise ValueError(f"Invalid host leaf name: {leaf.name!r}")
            observed_leaf_ids.add(int(match.group(1)))
        if observed_leaf_ids != expected_leaf_ids:
            raise ValueError(
                "Host leaves must be named consecutively H1,...,Hn. "
                f"Observed IDs: {sorted(observed_leaf_ids)}"
            )
        for leaf in leaves:
            leaf.name = leaf.name[1:]

        tree.name = str(num_leaves + 1)
        next_internal_id = num_leaves + 2
        for node in tree.traverse("preorder"):
            if node.is_root() or node.is_leaf():
                continue
            node.name = str(next_internal_id)
            next_internal_id += 1

        edges = []
        lengths = []
        # Preorder guarantees that the first edge starts at the real root.
        for node in tree.traverse("preorder"):
            if node.is_root():
                continue
            if node.dist is None or not np.isfinite(node.dist):
                raise ValueError(
                    f"Invalid branch length for node {node.name}: {node.dist}"
                )
            # Cophylo orders host speciation events by remaining branch length
            # and assumes an ultrametric tree, so branch lengths must be passed
            # through untouched.  Only guard against non-positive lengths, with
            # an epsilon small enough not to disturb the event ordering.
            edges.append(f"{node.up.name} {node.name}")
            lengths.append(f"{max(float(node.dist), 1e-9):.15g}")

        expected_edges = 2 * (num_leaves - 1)
        if len(edges) != expected_edges:
            raise ValueError(
                f"Expected {expected_edges} edges for a rooted binary tree with "
                f"{num_leaves} leaves, but generated {len(edges)}"
            )
        first_parent = int(edges[0].split()[0])
        if first_parent != num_leaves + 1:
            raise ValueError(
                f"First edge parent is {first_parent}; expected root "
                f"{num_leaves + 1}"
            )

        with open(edge_txt_path, "w") as edge_file:
            edge_file.write("\n".join(edges) + "\n")
        with open(length_txt_path, "w") as length_file:
            length_file.write("\n".join(lengths) + "\n")

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

    # -P: distribution of the number of hosts infected by each parasite lineage.
    # Entry i is the probability that a parasite infects exactly i+1 hosts, so
    # "1 0 0 ..." means strict specialists (no multiple associations).
    def save_specialist_distribution(self, tree, folder_path, tree_index): 
        k = len(tree.get_leaves())
        if self.host_distrib == "specialist":
            probs = np.zeros(k, dtype=np.float64)
            probs[0] = 1.0
        else:
            ranks = np.arange(1, k + 1, dtype=np.float64)
            probs = 1.0 / np.power(ranks, self.zipf_s)
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
    parser.add_argument("--regime", choices=sorted(REGIMES), default="medium",
                        help="Cospeciation/host-switch regime (sets -c and -s ranges)")
    parser.add_argument("--cosp_range", type=float, nargs=2, metavar=("MIN", "MAX"),
                        help="Override the regime's cospeciation probability range (-c)")
    parser.add_argument("--switch_range", type=float, nargs=2, metavar=("MIN", "MAX"),
                        help="Override the regime's host-switch rate range (-s)")
    parser.add_argument("--host_death", type=float,
                        help="Override the regime's host extinction rate (mu_H)")
    parser.add_argument("--sym_death", type=float,
                        help="Override the regime's symbiont extinction rate (mu_S, -m)")
    parser.add_argument("--min_parasites", type=int, default=5,
                        help="Minimum number of surviving parasite tips (cophylo -N)")
    parser.add_argument("--host_distrib", choices=["zipf", "specialist"], default="zipf",
                        help="Distribution of hosts per parasite (cophylo -P). "
                             "'specialist' forces one host per parasite.")
    parser.add_argument("--zipf_s", type=float, default=1.6,
                        help="Exponent of the zipf hosts-per-parasite distribution")
    parser.add_argument("--overwrite", action="store_true",
                        help="Delete previously generated folders in --output_dir first")
    args = parser.parse_args()

    regime = REGIMES[args.regime]
    cosp_range = tuple(args.cosp_range) if args.cosp_range else regime["cosp"]
    switch_range = tuple(args.switch_range) if args.switch_range else regime["switch"]
    host_death = args.host_death if args.host_death is not None else regime["host_death"]
    sym_death = args.sym_death if args.sym_death is not None else regime["sym_death"]

    if args.overwrite:
        clean_output_dirs(args.output_dir)

    total_start = time.time()

    # Generate Host Trees
    host_generator = GenerateHostTree(
        args.num_trees, args.min_leaves, args.max_leaves, args.output_dir,
        cosp_range=cosp_range, switch_range=switch_range,
        host_death=host_death, sym_death=sym_death,
        min_parasites=args.min_parasites,
        host_distrib=args.host_distrib, zipf_s=args.zipf_s,
    )
    print(f"Regime '{args.regime}': c ~ U{cosp_range}, s ~ U{switch_range}, "
          f"lambda_H={HOST_BIRTH_RATE}, mu_H={host_death}, "
          f"lambda_S={SYMBIONT_BIRTH_RATE}, mu_S={sym_death}, "
          f"T={host_generator.tmrca}, -N {args.min_parasites}, -P {args.host_distrib}"
          + (f" (k^-{args.zipf_s})" if args.host_distrib == "zipf" else ""))
    start = time.time()
    if not os.listdir(host_generator.host_tree_dir):
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
