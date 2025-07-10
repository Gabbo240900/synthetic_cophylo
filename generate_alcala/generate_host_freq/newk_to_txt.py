from ete3 import Tree

def newick_to_txt(newick_path, edge_txt_path, length_txt_path):
    tree = Tree(newick_path, format=1)
    edges = []
    lengths = []

    for i, node in enumerate(tree.traverse("postorder")):
        if not node.name:
            node.name = f"IH{i}"

    for node in tree.traverse("postorder"):
        if not node.is_root():
            parent = node.up
            edges.append(f"{parent.name} {node.name}")
            lengths.append(str(node.dist))

    with open(edge_txt_path, "w") as edge_file:
        edge_file.write("\n".join(edges))

    with open(length_txt_path, "w") as length_file:
        length_file.write("\n".join(lengths))

def txt_to_newick(edge_txt_path, length_txt_path, output_newick_path):
    from collections import defaultdict

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

    def attach_children(node_name):
        for child_name, _ in child_map.get(node_name, []):
            nodes[node_name].add_child(nodes[child_name])
            attach_children(child_name)

    attach_children(root_name)
    root.write(outfile=output_newick_path, format=5)

newick_to_txt("/Users/gabriele/Co-phyloformer/generate_host_freq/generated_trees/branches_lengths/host_tree_1.nwk", "host_edge_1.txt", "host_edge_lengths_1.txt")
txt_to_newick("host_edge_1.txt", "host_edge_lengths_1.txt", "host_tree_1.nwk")