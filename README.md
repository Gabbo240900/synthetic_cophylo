#  Similarities, differences and biases in cophylogenetic models for host-symbiont coevolution

This project contains the code to reproduce data shown in the paper ' Similarities, differences and biases in cophylogenetic models for host-symbiont coevolution'. The repository also contains a detailed jupter notebook with all the different anaylsis that have been presented in the experiment section of the aforementioned paper. 

# How to run the models 

Coala, treeducken, and Asymmetree do not require you to change input parameters to reproduce the data we show in our paper. 
For cophylo instead you need to run the generator three different times with the different input parameters. 

Below you find the code to run the respective models.

Make sure you have the reuiqred python packages from reuiqrments.txt and the required R libriareis in renv.lock
## Coala

python ..generate_coala/simulate_input_trees.py --num_trees 10000 	\
 --base_output_dir ./generated_trees \
 --num_threads 16 \
 --jar_path ./cophylogeny-ML/code/coala/TGLGenerator.jar

## Treeducken 

python generate_treeducken/simulate_input_files.py --num_trees 1000


## Asymmetree 

python generate_asymmetree/generate_asymmetree.py



## cophylo (Alcala's model)

Run the generator once per regime; `--regime` sets the cospeciation probability
(`-c`) and host-switch rate (`-s`) ranges, so no parameter has to be edited by hand.

    python generate_alcala/simulate_input_trees.py --num_trees 1000 \
    --min_leaves 15 --max_leaves 50 --regime medium \
    --output_dir ./generate_alcala/generated_trees_medium --overwrite

    python generate_alcala/simulate_input_trees.py --num_trees 1000 \
    --min_leaves 15 --max_leaves 50 --regime highCosp \
    --output_dir ./generate_alcala/generated_trees_highCosp --overwrite

    python generate_alcala/simulate_input_trees.py --num_trees 1000 \
    --min_leaves 15 --max_leaves 50 --regime highSwitch \
    --output_dir ./generate_alcala/generated_trees_highSwitch --overwrite

| regime | c (`-c`) | s (`-s`) | mu_H | mu_S (`-m`) |
| --- | --- | --- | --- | --- |
| `medium` | U(0.2, 0.4) | U(0.2, 0.4) | 0.45 | 0.6 |
| `highCosp` | U(0.7, 0.9) | U(0.0, 0.05) | 0.63 | 0.645 |
| `highSwitch` | U(0.05, 0.1) | U(0.5, 0.7) | 0.24 | 0.66 |

Shared across regimes: lambda_H = lambda_S = 0.7, T = TMRCA = 2 for both host and
symbiont, 15-50 extant host leaves, and a zipf hosts-per-symbiont distribution
f(k) proportional to k^-1.6.

Other options: `--cosp_range`/`--switch_range`/`--host_death`/`--sym_death` override
the regime, `--min_parasites`
sets cophylo's `-N` (minimum number of surviving parasite tips, default 5), and
`--host_distrib specialist` writes a `1 0 0 ...` hosts-per-parasite distribution
(`-P`) so that no parasite infects more than one host; the default `zipf`
distribution (exponent `--zipf_s`, default 1.6) allows generalists.

Note: cophylo stores its file names in 100-character buffers, so the generator
invokes it from inside `alcala_trees/` with short relative prefixes. Do not pass
absolute prefixes to `cophylo.out` by hand - it aborts without writing anything. 



