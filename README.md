# Cophylogeny simulators are not interchangeable: similarities, differences and structural biases in synthetic host–symbiont data

This repository contains the code and data to reproduce the results of the paper
'Cophylogeny simulators are not interchangeable: similarities, differences and structural biases in synthetic host–symbiont data'.
It includes:

- the generators of the four cophylogenetic models compared in the paper (Coala, treeducken,
  AsymmeTree and cophylo) and the synthetic datasets used in the paper;
- the real host–symbiont datasets used for comparison;
- the Jupyter notebooks with all the analyses presented in the experiments section of the paper.

All commands below are meant to be run from the root of the repository.

## Repository structure

| Path | Content |
| --- | --- |
| `generate_coala/` | Coala generator (`simulate_input_trees.py`) and the bundled `TGLGenerator.jar` |
| `generate_treeducken/` | treeducken generator (`simulate_input_files.py`, which calls `treeducken.r`) |
| `generate_asymmetree/` | AsymmeTree generator (`generate_asymmetree.py`) |
| `generate_alcala/` | cophylo (Alcala's model) generator (`simulate_input_trees.py`) and the cophylo C sources/binary in `software/` |
| `synthetic_datasets.py` | Location of the synthetic datasets of every generator and regime, used by the notebooks |
| `real_data/` | Real host–symbiont datasets (NEXUS) |
| `Analyze_trees.ipynb` | Tree and association-network comparison of the generators and the real data |
| `real_data_analysis/` | Tree-feature comparison of the real datasets with the synthetic ones |
| `regime_attainment/` | Distance of the generated event frequencies to the target regimes |
| `synthetic_data_analysis.ipynb` | Parsimonious reconciliation of every synthetic dataset with Capybara |
| `parsimonious_analysis/` | Analysis of the Capybara reconciliations |
| `requirements.txt` | Python dependencies |
| `renv.lock`, `renv/`, `.Rprofile` | R environment (managed with [renv](https://rstudio.github.io/renv/)) |

## Data

### Synthetic datasets

The repository contains the synthetic datasets used in the paper: about 1,000 datasets for each of
the four generators and each of the three regimes (**high cospeciation**, **high host switch** and
**medium**). Each dataset is a `.tgl` (NEXUS-like) file with the host tree, the symbiont tree, the
host–symbiont associations and the frequencies of the simulated events. They are stored in the
output folder of each generator:

| Generator | High cospeciation | High host switch | Medium |
| --- | --- | --- | --- |
| Coala | `generate_coala/generated_trees/high_cosp/Datasets/` | `.../high_switch/Datasets/` | `.../medium/Datasets/` |
| treeducken | `generate_treeducken/generated_trees_highCosp/Datasets/` | `.../generated_trees_highSwitch/Datasets/` | `.../generated_trees_medium/Datasets/` |
| AsymmeTree | `generate_asymmetree/generated_trees/high_cosp/Datasets/` | `.../high_switch/Datasets/` | `.../medium/Datasets/` |
| cophylo | `generate_alcala/generated_trees_highCosp/Dataset/` | `.../generated_trees_highSwitch/Dataset/` | `.../generated_trees_medium/Dataset/` |

The notebooks find these folders through `synthetic_datasets.py`. In the notebooks and in the
Capybara results the high-cospeciation regime of AsymmeTree is labelled `low_switch`.

`generate_treeducken/time_analysis/` and `generate_asymmetree/time_analysis/` contain additional
treeducken and AsymmeTree datasets simulated with a shorter (`small_time`) and a longer
(`big_time`) tree age.

### Real datasets

`real_data/` contains the real host–symbiont datasets in NEXUS format; `real_data/gene_species/`
contains four gene–species datasets.

## Setup

### Python

Python 3.11 is required. Create an environment and install the dependencies, e.g. with conda:

```bash
conda create -n cophylo python=3.11
conda activate cophylo
pip install -r requirements.txt
```

or with `venv`:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

This also installs [Capybara](https://github.com/Helio-Wang/Capybara-app)
(`capybara-cophylogeny`), used for the parsimonious reconciliations.

The steps below (R, Java and the cophylo binary) are needed only to generate new synthetic data;
the analysis notebooks need only Python.

### R (treeducken only)

R 4.5 is used. The R packages are pinned in `renv.lock`. `treeducken` is not on CRAN: it is
installed from our fork [Gabbo240900/treeducken](https://github.com/Gabbo240900/treeducken)
(pinned commit), which fixes the recording and counting of host-switch events
(`hs_mode = 'switch'`) of the original package. It is compiled from source, so a C++ compiler is
required (Xcode Command Line Tools on macOS, `build-essential` on Linux).

treeducken requests C++11, but the current RcppArmadillo headers need C++14 or newer. Tell R to
compile C++11 code with C++17 before restoring (only needed once, for the installation):

```bash
echo "CXX11STD = -std=gnu++17" > treeducken.Makevars
export R_MAKEVARS_USER="$PWD/treeducken.Makevars"
```

Then, from the same shell, start R from the repository root (the `.Rprofile` activates renv
automatically) and restore the environment:

```r
install.packages("renv")   # if renv is not installed yet
renv::restore()
```

`treeducken.r` activates the project's renv library by itself, so nothing else is needed when it is
called from the Python generator. `Rscript` must be on your `PATH`.

### Java (Coala only)

Coala's `TGLGenerator.jar` requires **Java 21 or newer**. Check with:

```bash
java -version
```

### cophylo binary (Alcala's model only)

The bundled binary `generate_alcala/software/bin/cophylo.out` is compiled for macOS on Apple Silicon (arm64).
On Linux or on Intel Macs recompile it (requires `gcc`):

```bash
gcc -O3 -o generate_alcala/software/bin/cophylo.out generate_alcala/software/sources/main_cophylo.c -lm
```

## Generating new data

The datasets of the paper are already in the repository; this section explains how to generate
new ones. **The generators write into the same folders as the paper's datasets**, so to keep them
either commit or copy the existing datasets first, or pass a different output folder (see below).

All generators use host trees of 15-50 extant leaves and a tree age (TMRCA) of 2. The commands
below reproduce the number of datasets used in the paper; change `--num_trees` to generate more
(or fewer) datasets.

### Coala

```bash
python generate_coala/simulate_input_trees.py --num_trees 10000 \
    --base_output_dir generate_coala/generated_trees \
    --num_threads 16
```

Output: `generate_coala/generated_trees/{high_cosp,high_switch,medium}/Datasets/`.
`--num_trees` is the number of host trees per regime: Coala does not produce a dataset for every
host tree (roughly 1 in 10 succeeds), so 10,000 host trees give about 1,000 datasets per regime.
The bundled `TGLGenerator.jar` is found automatically; use `--jar_path` (or the `COALA_TGL_JAR`
environment variable) to point to a different jar. `--num_threads` defaults to 75% of the CPUs.
Existing host trees are reused; pass `--regenerate` to create new ones.

### treeducken

```bash
python generate_treeducken/simulate_input_files.py --num_trees 1000
```

Output: `generate_treeducken/generated_trees_{highCosp,highSwitch,medium}/Datasets/`.

### AsymmeTree

```bash
python generate_asymmetree/generate_asymmetree.py --num_trees 1000
```

Output: `generate_asymmetree/generated_trees/{high_cosp,high_switch,medium}/Datasets/`
(change it with `--output_dir`).

### cophylo (Alcala's model)

Run the generator once per regime; `--regime` sets the cospeciation probability
(`-c`) and host-switch rate (`-s`) ranges, so no parameter has to be edited by hand.

```bash
python generate_alcala/simulate_input_trees.py --num_trees 1000 \
    --min_leaves 15 --max_leaves 50 --regime medium \
    --output_dir generate_alcala/generated_trees_medium --overwrite

python generate_alcala/simulate_input_trees.py --num_trees 1000 \
    --min_leaves 15 --max_leaves 50 --regime highCosp \
    --output_dir generate_alcala/generated_trees_highCosp --overwrite

python generate_alcala/simulate_input_trees.py --num_trees 1000 \
    --min_leaves 15 --max_leaves 50 --regime highSwitch \
    --output_dir generate_alcala/generated_trees_highSwitch --overwrite
```

Output: `generate_alcala/generated_trees_{medium,highCosp,highSwitch}/Dataset/`.
`--overwrite` deletes the previous datasets in `--output_dir`.

| regime | c (`-c`) | s (`-s`) | mu_H | mu_S (`-m`) |
| --- | --- | --- | --- | --- |
| `medium` | U(0.2, 0.4) | U(0.2, 0.4) | 0.45 | 0.6 |
| `highCosp` | U(0.7, 0.9) | U(0.0, 0.05) | 0.63 | 0.645 |
| `highSwitch` | U(0.05, 0.1) | U(0.5, 0.7) | 0.24 | 0.66 |

Shared across regimes: lambda_H = lambda_S = 0.7, T = TMRCA = 2 for both host and
symbiont, 15-50 extant host leaves, and a zipf hosts-per-symbiont distribution
f(k) proportional to k^-1.6.

Other options: `--cosp_range`/`--switch_range`/`--host_death`/`--sym_death` override
the regime, `--min_parasites` sets cophylo's `-N` (minimum number of surviving parasite
tips, default 5), and `--host_distrib specialist` writes a `1 0 0 ...` hosts-per-parasite
distribution (`-P`) so that no parasite infects more than one host; the default `zipf`
distribution (exponent `--zipf_s`, default 1.6) allows generalists.

Note: cophylo stores its file names in 100-character buffers, so the generator
invokes it from inside `alcala_trees/` with short relative prefixes. Do not pass
absolute prefixes to `cophylo.out` by hand - it aborts without writing anything.

### Using datasets in other folders

The notebooks read the datasets from the folders listed in `DATASET_DIRS` in
`synthetic_datasets.py`. To analyse datasets written elsewhere, edit the corresponding entries.

## Running the analyses

Start Jupyter with the Python environment activated:

```bash
jupyter lab
```

Each notebook finds the repository root by itself, so it can be opened from any folder. Tables and
figures are saved next to the notebook. Run the cells of each notebook in order.

| Notebook | What it does | Inputs |
| --- | --- | --- |
| `Analyze_trees.ipynb` | Tree statistics (leaves, branches, depth, balance indices), event counts and bipartite association-network metrics of the four generators, compared with the real datasets. Figures are saved in `comparison_trees/` and `comparison_networks/`. | synthetic and real datasets |
| `real_data_analysis/empirical_tree_feature_analysis.ipynb` | Standardised tree features (size, cherries, Sackin and Colless indices) of the real datasets and of a balanced sample of synthetic datasets; PCA and real-to-synthetic distances. | synthetic datasets, `real_data/*.nex` |
| `regime_attainment/target_region_graphics.ipynb` | Where the event frequencies of each generator fall with respect to the target regions R1 (high cospeciation), R2 (high host switch) and R3 (medium). | synthetic datasets |
| `synthetic_data_analysis.ipynb` | Runs Capybara (task 2, enumeration of the optimal reconciliations) on every synthetic dataset with the classical cost vectors <0,1,2,1> and <0,1,3,1> and with costs derived from the generating frequencies. Writes `synthetic_data_capybara_results.xlsx` and checkpoints in `capybara_checkpoints/`. | synthetic datasets |
| `parsimonious_analysis/capybara_downstream_analysis.ipynb` | Distance between the generating and the parsimonious event frequencies, and ambiguity of the solution space. | `synthetic_data_capybara_results.xlsx` |
| `parsimonious_analysis/capybara_optimal_profile_analysis.ipynb` | Event profiles of the optimal reconciliations, by generator and regime. | `synthetic_data_capybara_results.xlsx`, `synthetic_data_capybara_results_retry_asymmetree_300s.xlsx` |

The Capybara results of the paper are included in the repository root
(`synthetic_data_capybara_results.xlsx`, and `synthetic_data_capybara_results_retry_asymmetree_300s.xlsx`
with the AsymmeTree tasks rerun with a 300 s time limit), so the notebooks in
`parsimonious_analysis/` can be run without recomputing them.

To recompute them, run `synthetic_data_analysis.ipynb`. Processing all ~12,000 datasets takes
many hours: each Capybara task is stopped after `CAPYBARA_TIMEOUT_SECONDS` (60 s by default; the
time limit uses `SIGALRM`, so it works on macOS and Linux only). While it runs, the results
obtained so far are saved in `capybara_checkpoints/results_latest.csv`.
