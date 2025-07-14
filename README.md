# How to run the models 

## simulate with Alcala model high Switch
python simulate_input_trees.py --num_trees 10 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generate_host_freq/generated_trees/ 

(l = 0.7
m = 0.66
s = 0.7
c = 0.1
tmrca = 2
host_birth = 0.7
host_death = 0.24
)

## simulate with Alcala model high Cospeciation
python simulate_input_trees.py --num_trees 10 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generate_host_freq/generated_trees/ 

(l = 0.7
m = 0.615
s = 0.05
c = 1.0
tmrca = 2
host_birth = 0.7
host_death = 0.6
)


## simulate with Coala High switch - it was able to create only 3 trees
python simulate_input_trees.py --num_trees 100 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generated_trees/ \
--output_dir_tgl ./generated_trees/Datasets/ \
--jar_path ./cophylogeny-ML/code/coala/TGLGenerator.jar \
--num_threads 8

(pc = random(0.05,0.1)
ps = 0.7
pd = random.dirichlet * (1 - pc -ps) 
pl = 100-pc-ps-pd
host_birth = 0.7
host_death = 0.24)

## simulate with Coala High cospeciation 
python simulate_input_trees.py --num_trees 100 \
--min_leaves 15 --max_leaves 50 \
--output_dir ./generated_trees/ \
--output_dir_tgl ./generated_trees/Datasets/ \
--jar_path ./cophylogeny-ML/code/coala/TGLGenerator.jar \
--num_threads 8

(pc = random(0.7,1.0)
ps = 0.05
pd = random.dirichlet * (1 - pc -ps) 
pl = 100-pc-ps-pd
host_birth = 0.7
host_death = 0.63)



## Simulate with treeducken high_switch
python simulate_input_files.py \
  --h_lambda 0.7 0.7 \
  --c_lambda 0.1 0.1 \
  --s_lambda 0.7 0.7 \
  --s_her 1.4 1.4 \
  --num_trees 10

(h_mu = 0.24
s_mu = 0.66
t = 2)

## Simulate with treeducken high_cospeciation
python simulate_input_files.py \
  --h_lambda 0.7 0.7 \
  --c_lambda 1.4 1.4 \
  --s_lambda 0.7 0.7 \
  --s_her 0.05 0.05 \
  --num_trees 10

(h_mu = 0.63
s_mu = 0.645
t = 2)


# How the models work 

## ALCALA
From a given host tree, it generates a parasite tree. Works in C with a model that, in our examples, always fails to create associations because of memory allocation issues. It generates 5 files and takes as input 3 files that I then put together in a TGL file. 

The parameters for the parasite tree that I can input are: 
- l: Parasite speciation rate
- m: Parasite extinction rate
- c: Cospeciation ***probability***
- s: host switch rate
- t: tmrca (time to the most recent common ancestor)

The host parameters can be set when creating random host trees using the BioPython library with the Birth/Death model: 
- Host Birth Rate
- Host Death Rate
- Number of Leaves

The input files consist of:
- Host edges (.txt)
- Host edge lengths (.txt)
- Host Distributions (determines if parasite is a specialist or a generalist) -> list of num_leaves values that indicate the probability for any parasite node to associate with nth host leaf. (.txt)

The output files consist of:
- Parasite edges (.txt)
- Parasite edge lengths (.txt)
- Parasite number of nodes (leaves)(.txt)
- Parasite-Host Associations (.txt)
- Parasite Parameters (Cospeciation and Host Switch values) (.txt)

## COALA
From a given Host tree, it generates a Parasite tree with the relative associations and some coevolution parameters. It is written in Java and does not support dated trees (ignores host and parasite branch lengths).

The parameters I can set in the input are the following: 
- pc: Cospeciation probability 
- ps: Switch Probability
- pd: Duplication Probability
- Loss probability is measured as 100 - (pc+ps+pd)
- 
The host parameters can be set when creating random host trees using the BioPython library with the Birth/Death model: 
- Host Birth Rate
- Host Death Rate
- Number of Leaves

As input file for the Java generator, you need to give:
- Host tree in Newick format without branch lengths

The output of the TGL Generator:
- Dataset.tgl file (Nexus file) -> this file contains the host and parasite trees without branch lengths, the association mappings and some final statistics:
	- Number of leaves
	- Expected coevolution events (discrete)
	- Observed coevolution events (discrete)
	- Probabilities of the coevolution (input of Coala)
	- Frequencies (Observed distribution of coevolution events - continuous value)
	- Euclidean3
	- Euclidean4 

As a personal project, we attempted to implement a posterior branch length in Coala by utilizing the associations and branch lengths of the input Host tree. 
- We measure the depth of the tree and create a 'time variable' based on the depth (the older you are, the more time you have - the root has the smallest depth but the highest time).
- Based on this, we assign a time value to parasite nodes associated with a host node; it is the host value minus a tweak (more recent than the host).
- For internal nodes, we use a recursive function that computes the maximum time of the children and adds a small random increment (to make it older). 
- From the time assigned to parasite nodes, we do a reverse of the initial function (postorder traverse) to generate branch lengths. 

## TREEDUCKEN
This model is different from Coala and Alcala's project since it generates both the Host and Parasite tree without having to pass through BioPython. 
It is an R package (outdated and no anymore on CRAN). It has a problem with hs_mode that does not allow anymore to switch between only host_switches and only host_expansion. All the other parts work perfectly well. 

The parameters that we can set as input refer both to parasite and host tree: 
- Host speciation Rate
- Host Extinction Rate
- Parasite Speciation Rate
- Parasite Extinction Rate
- Cospeciation Rate
- Host Switch/Spread(expansion)

We do not have to give anything as input beside the parameters defined above 

The output of treeducken are some objects in R that we can decompose and obtain the following information: 
- Host and Parasite tree structure 
- Associations 
- Summary statistics: 
	- Cospeciation, 
	- Host Speciations,
	- Host Extinctions, 
	- Symbiont extinctions
	- Symbiont Speciations
	- Host Spread/Switches
	- Dispersal
	- Extirpations
	- Parafit Stat
	- Parafit P-value
All the statistics are discrete values and represent the number of times that the event happened during the coevolution (except for Dispersal and Extirpations, which are always N/A, and Parafit stat and P-value, which represent some metrics to understand how good the two trees coevolved). 

We can also set hs_mode which allows us to change between 'only switch', 'only expansion' or 'both'. In the latest version of the packet this works only with Both and only epxansion while only swithcing always returns N/A.

