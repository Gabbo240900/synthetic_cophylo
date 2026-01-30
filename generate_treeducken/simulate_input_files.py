import argparse
import subprocess
import os
import glob
import re 
import pandas as pd
import random
from ete3 import Tree
from tqdm import tqdm
import shutil

class GenerateTGLFiles:
    def __init__(self, h_lambda, h_mu, c_lambda, s_lambda, s_mu, s_her, time_to_sim, num_trees, base_out_dir, numbsim=1):
        self.h_lambda = h_lambda
        self.h_mu = h_mu
        self.c_lambda = c_lambda
        self.s_lambda = s_lambda
        self.s_mu = s_mu
        self.s_her = s_her
        self.time_to_sim = time_to_sim
        self.num_trees = num_trees
        self.numbsim = numbsim
        self.base_out_dir = base_out_dir

    def ensure_output_dirs_and_symlink(self):
        """
        Ensure the expected output subdirectories exist under base_out_dir and
        create/update a symlink named 'generated_trees' that points to base_out_dir.

        This is required because treeducken.r writes to 'generated_trees/...'.
        """
        # Create expected subdirectories
        for sub in ["host_tree", "symb_tree", "associations", "summaries", "Datasets"]:
            os.makedirs(os.path.join(self.base_out_dir, sub), exist_ok=True)

        # Create/update symlink in the script directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        link_path = os.path.join(script_dir, "generated_trees")

        # If link_path exists as a symlink, replace it
        if os.path.islink(link_path):
            os.unlink(link_path)
        elif os.path.exists(link_path):
            # If it's a real directory/file, do not delete it automatically (safety).
            raise RuntimeError(
                f"Path '{link_path}' exists and is not a symlink. "
                f"Please rename/remove it so the script can create a symlink "
                f"to '{self.base_out_dir}'."
            )

        os.symlink(os.path.abspath(self.base_out_dir), link_path)

    def run_r_script(self, sim_index):
        """Run the R script with parameters."""
        # treeducken.r writes into 'generated_trees/...', so ensure the symlink points
        # to this experiment's base_out_dir before running.
        self.ensure_output_dirs_and_symlink()

        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script_path = os.path.join(script_dir, "treeducken.r")

        args = [
            "Rscript",
            r_script_path,
            str(self.h_lambda),
            str(self.h_mu),
            str(self.c_lambda),
            str(self.s_lambda),
            str(self.s_mu),
            str(self.s_her),
            str(self.time_to_sim),
            str(sim_index)
        ]
        try:
            subprocess.run(
                args,
                check=True,
                cwd=script_dir,  # ensures symlink + relative paths resolve correctly
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError as e:
            print(f"Error running R script: {e}")
        except RuntimeError as e:
            print(f"Error preparing output paths: {e}")
    
    def generate_tgl_files(self, sim_index):
        # Find the most recent files for this simulation
        host_file = f"{self.base_out_dir}/host_tree/host_tree_{sim_index}.nwk"
        symb_file = f"{self.base_out_dir}/symb_tree/symb_tree_{sim_index}.nwk"
        assoc_file = f"{self.base_out_dir}/associations/association_{sim_index}.csv"
        summary_file = f"{self.base_out_dir}/summaries/summary_{sim_index}.csv"

        if not all(os.path.exists(f) for f in [host_file, symb_file, assoc_file, summary_file]):
            print(f"[WARN] One or more files missing for sim_index {sim_index}, skipping TGL generation.")
            return
        
        with open(host_file) as f:
            host_tree = f.read().strip()
        # Replace 'X' with 'H' in host_tree
        host_tree = host_tree.replace('X', 'H')
        with open(symb_file) as f:
            symb_tree = f.read().strip()
        # Replace 'X' with 'P' and 'S' with 'P' in symb_tree
        symb_tree = symb_tree.replace('X', 'P')
        symb_tree = symb_tree.replace('S', 'P')
        assoc_df = pd.read_csv(assoc_file)
        # Replace 'S' with 'P' in symbiont and host columns
        assoc_df['symbiont'] = assoc_df['symbiont'].str.replace('S', 'P')
        assoc_df['host'] = assoc_df['host'].str.replace('S', 'P')
        assoc_lines = [f'\t\t{row["symbiont"]}: {row["host"]}' for _, row in assoc_df.iterrows()]
        assoc_block = "\n".join(assoc_lines).rstrip(',')
        summary_df = pd.read_csv(summary_file, index_col=0).T.iloc[0]

        # Normalize event counts to proportions
        event_keys = [
            "Cospeciations",
            "Host_Spreads/switches",
            "Symbiont_Speciations",
            "Host_Extinctions",
            "Symbiont_Extinctions",
            'Host_Speciations'
        ]
        total_events = summary_df[event_keys].sum()
        for key in event_keys:
            summary_df[key] = summary_df[key] / total_events
            if summary_df[key] > 1.0:
                print(f"[WARNING] Normalized values exceed 1 for sim_index {sim_index}")

        host_leaf_count = len(Tree(host_tree, format=1).get_leaves())
        symb_leaf_count = len(Tree(symb_tree, format=1).get_leaves())

        summary_df["Host_num_leaves"] = host_leaf_count
        summary_df["Parasite_num_leaves"] = symb_leaf_count
        summary_df["Sim_time"] = self.time_to_sim
        summary_df['Total Events'] = total_events

        content = f"""#NEXUS
BEGIN HOST;
    TREE * Host1 = {host_tree};
ENDBLOCK;

BEGIN PARASITE;
    TREE * Para1 = {symb_tree};
ENDBLOCK;

BEGIN DISTRIBUTION;
    RANGE
{assoc_block}
    ;
END;

{summary_df.to_string(index=True, header=False)}
"""
        out_dir = os.path.join(self.base_out_dir, 'Datasets')
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f'Dataset{sim_index}.tgl')
        with open(out_path, 'w') as out_file:
            out_file.write(content)

def main():
    parser = argparse.ArgumentParser(description="Run Treeducken simulation and generate outputs.")
    parser.add_argument("--num_trees", type=int, required=True)

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    time_grid = [2]
    
    # INCREASE COSPECIATION VALUEs and switches etc... need to increase frequencies 

    experiments = [
        {
            "name": "highSwitch",
            "out_dir": os.path.join(script_dir, "time_analysis/small_time/high_switch"),
            "h_lambda": (0.7, 0.7),
            "c_lambda": (0.0, 0.05),
            "s_lambda": (0.7, 0.7),
            "s_her": (1.8, 2.2),
        },
        {
            "name": "highCosp",
            "out_dir": os.path.join(script_dir, "time_analysis/small_time/high_cosp"),
            "h_lambda": (0.7, 0.7),
            "c_lambda": (1.8, 2.2),
            "s_lambda": (0.7, 0.7),
            "s_her": (0.0, 0.05),
        },
        {
            "name": "medium",
            "out_dir": os.path.join(script_dir, "time_analysis/small_time/medium"),
            "h_lambda": (0.7, 0.7),
            "c_lambda": (0.5, 0.8),
            "s_lambda": (0.7, 0.7),
            "s_her": (0.5, 0.8),
        },
    ]

    for exp in experiments:
        print(f"Running experiment: {exp['name']}")
        os.makedirs(exp["out_dir"], exist_ok=True)

        sim_index = 1  # reset for each experiment

        for _ in tqdm(range(args.num_trees), desc=f"Simulating {exp['name']}"):
            h_lambda = random.uniform(*exp["h_lambda"])
            c_lambda = random.uniform(*exp["c_lambda"])
            s_lambda = random.uniform(*exp["s_lambda"])
            s_her = random.uniform(*exp["s_her"])

            h_mu = 0.3 * (h_lambda + c_lambda)
            s_mu = 0.3 * (s_lambda + c_lambda + s_her)

            for time_to_sim in time_grid:
                generator = GenerateTGLFiles(
                    h_lambda=h_lambda,
                    h_mu=h_mu,
                    c_lambda=c_lambda,
                    s_lambda=s_lambda,
                    s_mu=s_mu,
                    s_her=s_her,
                    time_to_sim=time_to_sim,
                    num_trees=args.num_trees,
                    base_out_dir=exp["out_dir"],
                    numbsim=1
                )
                generator.run_r_script(sim_index)
                generator.generate_tgl_files(sim_index)
                sim_index += 1

if __name__ == "__main__":
    main()
