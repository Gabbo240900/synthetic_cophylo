# How to run the models 


# Alcala

## simulate with Alcala model high Switch
python ..generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir generated_trees/ 


(l = 0.7
m = 0.66
s = random.uniform(0.5, 0.7)
c = random.uniform(0.05, 0.1)
tmrca = 2
host_birth = 0.7
host_death = 0.24
)

## simulate with Alcala model high Cospeciation
python ..generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generate_host_freq/generated_trees/ 

(l = 0.7
m = 0.645
s = random.uniform(0, 0.05)
c = random.uniform(0.7, 0.9)
tmrca = 2
host_birth = 0.7
host_death = 0.63
)

## simulate with Alcala model medium
python ..generate_alcala/simulate_input_trees.py --num_trees 1000 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generate_host_freq/generated_trees/ 

(l = 0.7
m = 0.6
s = random.uniform(0.2, 0.4)
c = random.uniform(0.2, 0.4)
tmrca = 2
host_birth = 0.7
host_death = 0.45
)



---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# COALA

python ..generate_coala/simulate_input_trees.py --num_trees 10000 	\
 --base_output_dir ./generated_trees \
 --num_threads 16 \
 --jar_path ./cophylogeny-ML/code/coala/TGLGenerator.jar


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Treeducken 

python ..generate_treeducken/simulate_input_files.py --num_trees 1000

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Asymmetree 

python ..generate_asymmetree/generate_asymmetree.py
