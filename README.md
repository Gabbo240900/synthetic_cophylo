#  Similarities, differences and biases in cophylogenetic models for host-symbiont coevolution

This project contains the code to reproduce data shown in the paper ' Similarities, differences and biases in cophylogenetic models for host-symbiont coevolution'. The repository also contains a detailed jupter notebook with all the different anaylsis that have been presented in the experiment section of the aforementioned paper. 

# How to run the models 

Coala, treeducken, and Asymmetree do not require you to change input parameters to reproduce the data we show in our paper. 
For cophylo you also run the generator three times, once per regime, selected with the `--regime` flag (no manual editing of the script is needed). 

Below you find the code to run the respective models.

Make sure you have the reuiqred python packages from reuiqrments.txt and the required R libriareis in renv.lock
## Coala

python ..generate_coala/simulate_input_trees.py --num_trees 10000 	\
 --base_output_dir ./generated_trees \
 --num_threads 16 \
 --jar_path ./cophylogeny-ML/code/coala/TGLGenerator.jar

## Treeducken 

python ..generate_treeducken/simulate_input_files.py --num_trees 1000


## Asymmetree 

python ..generate_asymmetree/generate_asymmetree.py



## cophylo (Alcala's model)

Run once per regime (`--regime` sets the ranges the cophylo `-c` cospeciation
probability and `-s` host-switch rate are drawn from for each tree):

python ../generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generated_trees --regime cospeciation

python ../generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generated_trees --regime mixed

python ../generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generated_trees --regime hostswitch

Available regimes are `cospeciation` (`-c` 0.8-1.0, `-s` 0.0-0.1), `mixed`
(`-c` 0.2-0.4, `-s` 0.2-0.4), `hostswitch` (`-c` 0.0-0.2, `-s` 0.8-1.0) and
`random` (both drawn from cophylo's uniform priors). `--cospec_range` and
`--switch_range` override the regime.

Other options:

* `--min_parasite_tips` (cophylo `-N`, default 5): minimum number of surviving
  parasite tips demanded of each simulation.
* `--host_dist` (cophylo `-P`, default `zipf`): distribution of the number of
  hosts per new parasite lineage. `zipf` leaves roughly half the probability
  mass on two or more hosts; `specialist` writes `1 0 0 0`, so each new lineage
  starts on a single host and multi-host associations arise only from host
  switches.

The script must be run from the `generate_alcala` directory, since it invokes
`software/bin/cophylo.out` by relative path. 



