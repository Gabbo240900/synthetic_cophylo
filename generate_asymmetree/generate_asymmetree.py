import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick
from tralda.datastructures import Tree, LCA
from collections import Counter
import re
import os
import csv 
import random

base_path = "/Users/gabriele/synthetic_cophylo/generate_asymmetree/generated_trees"
num_trees = 1000
species_pattern = re.compile(r'(?<=\(|,)(\d+):')
gene_pattern = re.compile(r'(?<=[\(|,|)])(\d+)(?=(?:<|:))')
assoc_pattern = re.compile(r'<(\d+(?:-\d+)*)>')


event_name_map = {
    "L": "Loss",
    "S": "Speciation",
    "H": "Horizontal Gene Transfer",
    "D": "Duplication"
}


# High cosp 
for i in range(num_trees):
    os.makedirs(f"{base_path}/high_cosp/species_trees", exist_ok=True)
    os.makedirs(f"{base_path}/high_cosp/gene_trees", exist_ok=True)
    os.makedirs(f"{base_path}/high_cosp/associations", exist_ok=True)
    os.makedirs(f"{base_path}/high_cosp/reconciliations", exist_ok=True)
    os.makedirs(f"{base_path}/high_cosp/Datasets", exist_ok=True)
    num_leaves = random.randint(15, 50)
    s = te.species_tree_n_age(n=num_leaves, model='BDP', age=2, birth_rate=0.7, death_rate=0.63)
    #DECIDE HOW WE CHOOSE DUPLICATION AND LOSS RATES. 
    hgt_rate = random.uniform(0, 0.05)
    g = te.dated_gene_tree(s, dupl_rate=0.5, loss_rate=0.5, hgt_rate=hgt_rate) 
    s_nwk = to_newick(s)
    g_nwk = to_newick(g)
    s_nwk = re.sub(r'(?<=[(,)])(?!(H))(\d+):', r'H\2:', s_nwk)
    g_nwk = gene_pattern.sub(lambda m: f"P{m.group(1)}", g_nwk)
    g_nwk = assoc_pattern.sub(lambda m: "<" + re.sub(r'\d+', lambda n: "H" + n.group(0), m.group(1)) + ">", g_nwk)
    s_path = f"{base_path}/high_cosp/species_trees/species_tree_{i+1}.nwk"
    g_path = f"{base_path}/high_cosp/gene_trees/gene_tree_{i+1}.nwk"
    with open(s_path, "w") as f:
        f.write(s_nwk + "\n")
        
    associations = []
    for match in re.finditer(r'(P\d+)<([^>]+)>', g_nwk):
        p_node = match.group(1)
        h_nodes = match.group(2).split('-')
        for h in h_nodes:
            associations.append((p_node, h))

    # Save associations CSV
    assoc_path = f"{base_path}/high_cosp/associations/associations_{i+1}.csv"
    with open(assoc_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["parasite_leaf:", "host_leaf"])
        writer.writerows(associations)

    # Remove all angle brackets and their contents from g_nwk before writing
    g_nwk_clean = re.sub(r'<[^>]*>', '', g_nwk)
    with open(g_path, "w") as f:
        f.write(g_nwk_clean + "\n")
        
        
    event_counts = Counter(
        getattr(node, "event", None)
        for node in g.postorder()
        if getattr(node, "event", None)
    )
    gene_leaf_count = sum(1 for _ in g.leaves())
    speciation_events = event_counts.get("S", 0)
    adjusted_speciation = speciation_events - gene_leaf_count
    rec_path = f"{base_path}/high_cosp/reconciliations/reconciliation_{i+1}.txt"
    with open(rec_path, "w") as f:
        for event, count in event_counts.items():
            long_name = event_name_map.get(event, event)
            if event == "S":
                f.write(f"Speciation: {adjusted_speciation}\n")
            else:
                f.write(f"{long_name}: {count}\n")
    
    dataset_path = f"{base_path}/high_cosp/Datasets/Dataset{i+1}.tgl"
    with open(dataset_path, "w") as f:
        f.write("#NEXUS\n")

        # Host tree block
        f.write("BEGIN HOST;\n")
        f.write(f"\tTREE * Host1 = {s_nwk}\n")
        f.write("ENDBLOCK;\n\n")

        # Symbiont tree block
        f.write("BEGIN PARASITE;\n")
        f.write(f"\tTREE * Para1 = {g_nwk_clean}\n")
        f.write("ENDBLOCK;\n\n")

        # Associations block
        f.write("BEGIN DISTRIBUTION;\n")
        f.write('\tRANGE\n')
        for p_node, h_node in associations:
            f.write(f"\t\t{p_node}: {h_node}\n")
        f.write("END;\n\n")

        # Reconciliation block
        for event, count in event_counts.items():
            long_name = event_name_map.get(event, event)
            if event == "S":
                f.write(f"Speciation: {adjusted_speciation}\n")
            else:
                f.write(f"{long_name}: {count}\n")
        
print('Done generating high cosp species and gene trees.')

# High Switch 
for i in range(num_trees):
    os.makedirs(f"{base_path}/high_switch/species_trees", exist_ok=True)
    os.makedirs(f"{base_path}/high_switch/gene_trees", exist_ok=True)
    os.makedirs(f"{base_path}/high_switch/associations", exist_ok=True)
    os.makedirs(f"{base_path}/high_switch/reconciliations", exist_ok=True)
    os.makedirs(f"{base_path}/high_switch/Datasets", exist_ok=True)
    num_leaves = random.randint(15, 50)
    s = te.species_tree_n_age(n=num_leaves,model='BDP', age=2, birth_rate=0.7, death_rate=0.24)
    #DECIDE HOW WE CHOOSE DUPLICATION AND LOSS RATES. 
    hgt_rate = random.uniform(1.4, 1.7)
    g = te.dated_gene_tree(s, dupl_rate=0.5, loss_rate=0.5, hgt_rate=hgt_rate) 
    s_nwk = to_newick(s)
    g_nwk = to_newick(g)
    s_nwk = re.sub(r'(?<=[(,)])(?!(H))(\d+):', r'H\2:', s_nwk)
    g_nwk = gene_pattern.sub(lambda m: f"P{m.group(1)}", g_nwk)
    g_nwk = assoc_pattern.sub(lambda m: "<" + re.sub(r'\d+', lambda n: "H" + n.group(0), m.group(1)) + ">", g_nwk)
    s_path = f"{base_path}/high_switch/species_trees/species_tree_{i+1}.nwk"
    g_path = f"{base_path}/high_switch/gene_trees/gene_tree_{i+1}.nwk"
    with open(s_path, "w") as f:
        f.write(s_nwk + "\n")
        
    associations = []
    for match in re.finditer(r'(P\d+)<([^>]+)>', g_nwk):
        p_node = match.group(1)
        h_nodes = match.group(2).split('-')
        for h in h_nodes:
            associations.append((p_node, h))

    # Save associations CSV
    assoc_path = f"{base_path}/high_switch/associations/associations_{i+1}.csv"
    with open(assoc_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["parasite_leaf:", "host_leaf"])
        writer.writerows(associations)

    # Remove all angle brackets and their contents from g_nwk before writing
    g_nwk_clean = re.sub(r'<[^>]*>', '', g_nwk)
    with open(g_path, "w") as f:
        f.write(g_nwk_clean + "\n")
        
        
    event_counts = Counter(
        getattr(node, "event", None)
        for node in g.postorder()
        if getattr(node, "event", None)
    )
    gene_leaf_count = sum(1 for _ in g.leaves())
    speciation_events = event_counts.get("S", 0)
    adjusted_speciation = speciation_events - gene_leaf_count
    rec_path = f"{base_path}/high_switch/reconciliations/reconciliation_{i+1}.txt"
    with open(rec_path, "w") as f:
        for event, count in event_counts.items():
            long_name = event_name_map.get(event, event)
            if event == "S":
                f.write(f"Speciation: {adjusted_speciation}\n")
            else:
                f.write(f"{long_name}: {count}\n")
    
    dataset_path = f"{base_path}/high_switch/Datasets/Dataset{i+1}.tgl"
    with open(dataset_path, "w") as f:
        f.write("#NEXUS\n")

        # Host tree block
        f.write("BEGIN HOST;\n")
        f.write(f"\tTREE * Host1 = {s_nwk}\n")
        f.write("ENDBLOCK;\n\n")

        # Symbiont tree block
        f.write("BEGIN PARASITE;\n")
        f.write(f"\tTREE * Para1 = {g_nwk_clean}\n")
        f.write("ENDBLOCK;\n\n")

        # Associations block
        f.write("BEGIN DISTRIBUTION;\n")
        f.write('\tRANGE\n')
        for p_node, h_node in associations:
            f.write(f"\t\t{p_node}: {h_node}\n")
        f.write("END;\n\n")

        # Reconciliation block
        for event, count in event_counts.items():
            long_name = event_name_map.get(event, event)
            if event == "S":
                f.write(f"Speciation: {adjusted_speciation}\n")
            else:
                f.write(f"{long_name}: {count}\n")
        
print('Done generating high switch species and gene trees.')

# Medium
for i in range(num_trees):
    os.makedirs(f"{base_path}/medium/species_trees", exist_ok=True)
    os.makedirs(f"{base_path}/medium/gene_trees", exist_ok=True)
    os.makedirs(f"{base_path}/medium/associations", exist_ok=True)
    os.makedirs(f"{base_path}/medium/reconciliations", exist_ok=True)
    os.makedirs(f"{base_path}/medium/Datasets", exist_ok=True)
    num_leaves = random.randint(15, 50)
    s = te.species_tree_n_age(n=num_leaves,model='BDP', age=2, birth_rate=0.7, death_rate=0.45)
    #DECIDE HOW WE CHOOSE DUPLICATION AND LOSS RATES. 
    hgt_rate = random.uniform(0.5, 0.8)
    g = te.dated_gene_tree(s, dupl_rate=0.5, loss_rate=0.5, hgt_rate=hgt_rate) 
    s_nwk = to_newick(s)
    g_nwk = to_newick(g)
    s_nwk = re.sub(r'(?<=[(,)])(?!(H))(\d+):', r'H\2:', s_nwk)
    g_nwk = gene_pattern.sub(lambda m: f"P{m.group(1)}", g_nwk)
    g_nwk = assoc_pattern.sub(lambda m: "<" + re.sub(r'\d+', lambda n: "H" + n.group(0), m.group(1)) + ">", g_nwk)
    s_path = f"{base_path}/medium/species_trees/species_tree_{i+1}.nwk"
    g_path = f"{base_path}/medium/gene_trees/gene_tree_{i+1}.nwk"
    with open(s_path, "w") as f:
        f.write(s_nwk + "\n")
        
    associations = []
    for match in re.finditer(r'(P\d+)<([^>]+)>', g_nwk):
        p_node = match.group(1)
        h_nodes = match.group(2).split('-')
        for h in h_nodes:
            associations.append((p_node, h))

    # Save associations CSV
    assoc_path = f"{base_path}/medium/associations/associations_{i+1}.csv"
    with open(assoc_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["parasite_leaf:", "host_leaf"])
        writer.writerows(associations)

    # Remove all angle brackets and their contents from g_nwk before writing
    g_nwk_clean = re.sub(r'<[^>]*>', '', g_nwk)
    with open(g_path, "w") as f:
        f.write(g_nwk_clean + "\n")
        
        
    event_counts = Counter(
        getattr(node, "event", None)
        for node in g.postorder()
        if getattr(node, "event", None)
    )
    gene_leaf_count = sum(1 for _ in g.leaves())
    speciation_events = event_counts.get("S", 0)
    adjusted_speciation = speciation_events - gene_leaf_count
    rec_path = f"{base_path}/medium/reconciliations/reconciliation_{i+1}.txt"
    with open(rec_path, "w") as f:
        for event, count in event_counts.items():
            long_name = event_name_map.get(event, event)
            if event == "S":
                f.write(f"Speciation: {adjusted_speciation}\n")
            else:
                f.write(f"{long_name}: {count}\n")
    
    dataset_path = f"{base_path}/medium/Datasets/Dataset{i+1}.tgl"
    with open(dataset_path, "w") as f:
        f.write("#NEXUS\n")

        # Host tree block
        f.write("BEGIN HOST;\n")
        f.write(f"\tTREE * Host1 = {s_nwk}\n")
        f.write("ENDBLOCK;\n\n")

        # Symbiont tree block
        f.write("BEGIN PARASITE;\n")
        f.write(f"\tTREE * Para1 = {g_nwk_clean}\n")
        f.write("ENDBLOCK;\n\n")

        # Associations block
        f.write("BEGIN DISTRIBUTION;\n")
        f.write('\tRANGE\n')
        for p_node, h_node in associations:
            f.write(f"\t\t{p_node}: {h_node}\n")
        f.write("END;\n\n")

        # Reconciliation block
        for event, count in event_counts.items():
            long_name = event_name_map.get(event, event)
            if event == "S":
                f.write(f"Speciation: {adjusted_speciation}\n")
            else:
                f.write(f"{long_name}: {count}\n")
        
print('Done generating medium species and gene trees.')