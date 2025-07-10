### Simulate with treeducken

python simulate_input_files.py \
  --h_lambda 0.5 1.2 \
  --c_lambda 0.7 1.6 \
  --s_lambda 0.5 1.2 \
  --s_her 0.05 0.3 \
  --num_trees 200




### simulate with Coala
python simulate_input_trees.py --num_trees 500 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generated_trees/ \
--output_dir_tgl ./generated_trees/Datasets/ \
--jar_path ./cophylogeny-ML/code/coala/TGLGenerator.jar \
--num_threads 8

### simulate with Alcala model 
python simulate_input_trees.py --num_trees 100 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generate_host_freq/generated_trees/ 
