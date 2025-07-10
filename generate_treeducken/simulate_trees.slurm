#!/bin/bash

#SBATCH --job-name=simulate_trees

#SBATCH --partition=prepost 

#SBATCH --output=simulate_trees%j.out 

#SBATCH --error=simulate_trees%j.err  

#SBATCH --nodes=1

#SBATCH --ntasks-per-node=48 

#SBATCH --time=20:00:00 

# WORKING DIRECTORY 
cd $ALL_CCFRWORK/Co-phyloformer/generate_treeducken/


# ENVIRONMENT SETUP

module purge

module load r/4.4.1

# Activate R environment for treeducken
export R_LIBS_USER=$ALL_CCFRWORK/r_envs/treeducken_env/lib

export MPLCONFIGDIR=$ALL_CCFRWORK/tmp/mpl_cache
export TRANSFORMERS_CACHE=$ALL_CCFRWORK/tmp/transformers_cache

module load miniforge/24.9.0

conda activate myenv 


# ECHO

set -x 

# STEP 0: remove former datasets
rm -rf generated_trees/associations/*
rm -rf generated_trees/Datasets/*
rm -rf generated_trees/host_trees/*
rm -rf generated_trees/summaries/*
rm -rf generated_trees/symb_tree/*

# STEP 1: generate host trees
python simulate_input_files.py \
  --h_lambda 0.5 1.2 \
  --c_lambda 0.7 1.6 \
  --s_lambda 0.5 1.2 \
  --s_her 0.05 0.3 \
  --num_trees 20

python alisim.py generated_trees/Datasets \
--substitution LG \
--gamma GC \
--iqtree ../bin/bin_linux/iqtree_2.2.0 \
--length 500 \
--max-attempts 1 \
--allow-duplicate-sequences 





