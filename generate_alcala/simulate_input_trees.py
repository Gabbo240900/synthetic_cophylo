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
# Trees where cophylo could not honour -N, or reported a broken host association
unmet_nmin = []

# Cophylogenetic regimes: (range for cophylo -c, range for cophylo -s).
# -c is the cospeciation probability given a host speciation, -s the host-switch
# rate. The cospeciation/host-switch endpoints mirror the two regimes shipped in
# the cophylo software's own run_examples.sh; 'random' reproduces the uniform
# priors cophylo uses when -c and -s are left unset.
REGIMES = {
    "cospeciation": ((0.8, 1.0), (0.0, 0.1)),
    "mixed":        ((0.2, 0.4), (0.2, 0.4)),
    "hostswitch":   ((0.0, 0.2), (0.8, 1.0)),
    "random":       ((0.0, 1.0), (0.0, 1.0)),
}


# Clean existing folders if present
for folder in ["alcala_trees", "Dataset", "host_trees", "parasite_trees"]:
    full_path = os.path.join("generated_trees", folder)
    if os.path.exists(full_path):
        shutil.rmtree(full_path)

class GenerateHostTree:
    """Class to generate host trees and their corresponding event frequency files."""
    
    def __init__(self, num_trees, min_leaves, max_leaves, output_dir, tmrca=2,
                 cospec_range=(0.2, 0.4), switch_range=(0.2, 0.4),
                 min_parasite_tips=5, host_dist="zipf", zipf_exponent=1.6):
        self.num_trees = num_trees  
        self.min_leaves = min_leaves  
        self.max_leaves = max_leaves  
        self.output_dir = output_dir  
        self.tmrca = tmrca
        # Simulation regime: -c (cospeciation probability) and -s (host-switch rate).
        self.cospec_range = cospec_range
        self.switch_range = switch_range
        # -N: minimum number of surviving parasite tips demanded of cophylo.
        self.min_parasite_tips = min_parasite_tips
        # -P: distribution of the number of hosts per new parasite lineage.
        self.host_dist = host_dist
        self.zipf_exponent = zipf_exponent

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
        s = random.uniform(*self.switch_range)
        c = random.uniform(*self.cospec_range)

        cmd = [
            "software/bin/cophylo.out",
            "-l", "0.7",
            "-m", "0.6",
            '-c', str(c),
            '-s', str(s),
            "-t", str(self.tmrca),
            "-i", input_prefix,
            "-o", output_prefix,
            "-N", str(self.min_parasite_tips),
            "-S", str(random.randint(1, 1000000)),
            "-P", distrib_path
        ]
        try:
            proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
            # cophylo reports an unmet -N (and any surviving association problem)
            # on stderr; without this the dataset would be written anyway and the
            # shortfall would go unnoticed.
            if proc.stderr.strip():
                for line in proc.stderr.strip().splitlines():
                    print(f"[cophylo tree {host_tree_index}] {line}")
                if "Could not reach Nmin" in proc.stderr or "invalid host edge index" in proc.stderr:
                    unmet_nmin.append(host_tree_index)
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


    @staticmethod
    def unplant(tree):
        """Collapse a leading chain of single-child nodes and return the real root.

        AsymmeTree emits planted trees (root -> single child -> ...). cophylo's
        edge-list format has no room for a stem edge: it derives the number of
        hosts from the first node id and then reads exactly 2*(n-1) edges, so a
        stem would push one real edge past the end of what it reads.
        """
        while len(tree.children) == 1:
            child = tree.children[0]
            child.dist = 0.0
            child.detach()
            tree = child
        return tree

    def generate_random_tree(self, num_leaves, prefix="H"):
        s = te.species_tree_n_age(n=num_leaves, model='BDP', age=self.tmrca, birth_rate=0.7, death_rate=0.45)
        s_nwk = to_newick(s). strip()
        if not s_nwk.endswith(";"):
            s_nwk += ";"
        s_nwk = re.sub(r'(?<=[(,)])(?!(H))(\d+):', r'H\2:', s_nwk)
        tree = Tree(s_nwk, format=1)
        # AsymmeTree returns a *planted* tree: the root carries a single child
        # (the stem edge). cophylo expects an unplanted rooted binary tree with
        # 2n-1 nodes and 2n-2 edges, so drop the stem before doing anything else.
        tree = self.unplant(tree)
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
        """Write the host tree as a cophylo edge list, in the R ``phylo`` convention.

        cophylo reads the first integer of ``host_edge.txt`` and takes it to be
        the root, deriving the number of hosts as ``root_id - 1`` and then
        reading exactly ``2*(n-1)`` edges (see ``main_cophylo.c``: ``nh`` is read
        from the first token, and ``H_edge[0]`` is matched against edge parents
        to find the two branches descending from the root).

        The tree must therefore be numbered leaves ``1..n``, root ``n+1``,
        remaining internal nodes ``n+2..2n-1``, and the edges must be written so
        that the first one starts at the root. Writing edges in postorder puts a
        cherry on the first line instead, which makes cophylo mistake that cherry
        for the whole host tree.

        Leaf ids are taken from the ``H<k>`` names rather than from traversal
        order, because ``create_tgl_dataset`` maps a host id ``k`` reported by
        cophylo back to the leaf named ``H<k>``.
        """
        tree = Tree(newick_path, format=5)
        tree = self.unplant(tree)
        leaves = list(tree.iter_leaves())
        num_leaves = len(leaves)

        if num_leaves < 2:
            raise ValueError(f"{newick_path}: need at least 2 host leaves, got {num_leaves}")

        # cophylo assumes a rooted binary host tree.
        for node in tree.traverse():
            if not node.is_leaf() and len(node.children) != 2:
                raise ValueError(
                    f"{newick_path}: host tree is not binary: node {node.name!r} "
                    f"has {len(node.children)} children"
                )

        # Leaves must be named H1..Hn so that ids 1..n are unambiguous.
        observed = set()
        for leaf in leaves:
            match = re.fullmatch(r"H(\d+)", leaf.name or "")
            if not match:
                raise ValueError(f"{newick_path}: invalid host leaf name {leaf.name!r}")
            observed.add(int(match.group(1)))
        if observed != set(range(1, num_leaves + 1)):
            raise ValueError(
                f"{newick_path}: host leaves must be named H1..H{num_leaves}; "
                f"observed ids {sorted(observed)}"
            )

        # Strip the H prefix: cophylo wants bare integer node ids.
        for leaf in leaves:
            leaf.name = leaf.name[1:]

        # Root is n+1; the other internal nodes follow it.
        tree.name = str(num_leaves + 1)
        next_internal_id = num_leaves + 2
        for node in tree.traverse("preorder"):
            if node.is_root() or node.is_leaf():
                continue
            node.name = str(next_internal_id)
            next_internal_id += 1

        edges = []
        lengths = []
        # Preorder guarantees the first edge descends from the root.
        for node in tree.traverse("preorder"):
            if node.is_root():
                continue
            dist = node.dist
            if dist is None or not np.isfinite(dist) or dist > 1e6:
                raise ValueError(
                    f"{newick_path}: invalid branch length {dist!r} on node {node.name}"
                )
            # cophylo needs strictly positive branch lengths; keep edges and
            # lengths strictly 1:1 so line i of one file matches line i of the other.
            edges.append(f"{node.up.name} {node.name}")
            lengths.append(f"{max(float(dist), 0.01):.15g}")

        expected_edges = 2 * (num_leaves - 1)
        if len(edges) != expected_edges:
            raise ValueError(
                f"{newick_path}: expected {expected_edges} edges for a rooted binary "
                f"tree on {num_leaves} leaves, wrote {len(edges)}"
            )
        first_parent = int(edges[0].split()[0])
        if first_parent != num_leaves + 1:
            raise ValueError(
                f"{newick_path}: first edge starts at {first_parent}, expected root "
                f"{num_leaves + 1}; cophylo would infer the wrong root and host count"
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

    def save_specialist_distribution(self, tree, folder_path, tree_index):
        """Write the -P file: P(a new parasite lineage infects exactly k hosts), k=1,2,...

        ``host_dist="specialist"`` writes ``1 0 0 0``, i.e. every new lineage
        infects a single host and no multi-host associations are produced.
        ``host_dist="zipf"`` writes a Zipf distribution over 1..n_hosts, which
        leaves roughly half the mass on k>=2 (mean ~3 hosts per lineage), so it
        is a generalist-leaning prior despite the historical function name.
        """
        if self.host_dist == "specialist":
            values = ["1", "0", "0", "0"]
        else:
            k = len(tree.get_leaves())
            ranks = np.arange(1, k + 1, dtype=np.float64)
            probs = 1.0 / np.power(ranks, self.zipf_exponent)
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
    parser.add_argument(
        "--regime", type=str, default="mixed", choices=sorted(REGIMES),
        help="Cophylogenetic regime, i.e. the ranges the cophylo -c (cospeciation "
             "probability) and -s (host-switch rate) are drawn from per tree. "
             + "; ".join(f"{k}: c={v[0]}, s={v[1]}" for k, v in sorted(REGIMES.items()))
    )
    parser.add_argument("--cospec_range", type=float, nargs=2, metavar=("MIN", "MAX"),
                        help="Override the regime's -c range")
    parser.add_argument("--switch_range", type=float, nargs=2, metavar=("MIN", "MAX"),
                        help="Override the regime's -s range")
    parser.add_argument("--min_parasite_tips", type=int, default=5,
                        help="cophylo -N: minimum number of surviving parasite tips")
    parser.add_argument("--host_dist", type=str, default="zipf", choices=["zipf", "specialist"],
                        help="cophylo -P: distribution of hosts per new parasite lineage. "
                             "'specialist' writes '1 0 0 0' (no multi-host associations); "
                             "'zipf' keeps the generalist-leaning Zipf prior")
    parser.add_argument("--zipf_exponent", type=float, default=1.6,
                        help="Exponent of the Zipf host-per-parasite distribution")
    args = parser.parse_args()

    cospec_range = tuple(args.cospec_range) if args.cospec_range else REGIMES[args.regime][0]
    switch_range = tuple(args.switch_range) if args.switch_range else REGIMES[args.regime][1]
    if args.min_parasite_tips < 2:
        parser.error("--min_parasite_tips must be at least 2: a parasite tree with a "
                     "single tip carries no cophylogenetic signal")
    print(f"Regime '{args.regime}': -c in {cospec_range}, -s in {switch_range}, "
          f"-N {args.min_parasite_tips}, -P {args.host_dist}")

    total_start = time.time()

    # Generate Host Trees
    host_generator = GenerateHostTree(
        args.num_trees, args.min_leaves, args.max_leaves, args.output_dir,
        cospec_range=cospec_range, switch_range=switch_range,
        min_parasite_tips=args.min_parasite_tips,
        host_dist=args.host_dist, zipf_exponent=args.zipf_exponent,
    )
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

    if unmet_nmin:
        print(f"\ncophylo reported warnings for {len(unmet_nmin)} tree(s) "
              f"(unmet -N and/or broken host associations): {unmet_nmin}")