# How to run the models 

Coala, treeducken, and Asymmetree do not require you to change input parameters to reproduce the data we show in our paper. 
For cophylo instead you need to run the generator three different times with the different input parameters. 

Below you find the code to run the respective models.

Make sure you have the reuiqred python paclages from reuiqrments.txt and the required R libriareis in renv.lock
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

python ..generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generate_host_freq/generated_trees/ 



