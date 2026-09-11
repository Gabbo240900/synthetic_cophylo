# Cophylogeny simulators are not interchangeable: similarities, differences and structural biases in synthetic host–symbiont data

This repository contains the code to reproduce the data shown in the paper
'Cophylogeny simulators are not interchangeable: similarities, differences and structural biases in synthetic host–symbiont data'.
It includes the generators for the four cophylogenetic models compared in the paper
(Coala, treeducken, AsymmeTree and cophylo) and a Jupyter notebook,
`Analyze_trees.ipynb`, with all the analyses presented in the experiments section of the paper.

All commands below are meant to be run from the root of the repository.

## Repository structure

| Path | Content |
| --- | --- |
| `generate_coala/` | Coala generator (`simulate_input_trees.py`) and the bundled `TGLGenerator.jar` |
| `generate_treeducken/` | treeducken generator (`simulate_input_files.py`, which calls `treeducken.r`) |
| `generate_asymmetree/` | AsymmeTree generator (`generate_asymmetree.py`) |
| `generate_alcala/` | cophylo (Alcala's model) generator (`simulate_input_trees.py`) and the cophylo C sources/binary in `software/` |
| `real_data/` | Real host-symbiont datasets (NEXUS) used for comparison |
| `Analyze_trees.ipynb` | Analysis notebook |
| `requirements.txt` | Python dependencies |
| `renv.lock`, `renv/`, `.Rprofile` | R environment (managed with [renv](https://rstudio.github.io/renv/)) |

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

## Generating the data

Each generator produces the three regimes used in the paper: **high cospeciation**, **high host switch**
and **medium**. All generators use host trees of 15-50 extant leaves and a tree age (TMRCA) of 2.
The commands below reproduce the number of datasets used in the paper; change `--num_trees` to
generate more (or fewer) datasets. Every dataset is written as a `.tgl` (NEXUS-like) file containing
the host tree, the symbiont tree and the host-symbiont associations.

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

## Running the analysis notebook

1. Generate the data of all four models as described above, keeping the default output folders
   (the notebook reads the datasets from those paths).
2. Start Jupyter **from the repository root**, with the Python environment activated:

   ```bash
   jupyter lab Analyze_trees.ipynb
   ```

3. Run the cells in order. The first code cell checks that the working directory is the repository
   root and creates the output folders. Figures are saved in `comparison_trees/` and
   `comparison_networks/`.
