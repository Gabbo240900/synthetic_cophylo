###############################################################################
#### Script to run some example simulations of host-parasite cophylogenies ####
###############################################################################

#Small example with cospeciation
bin/cophylo.out -l 0.7 -m 0.615 -t 2 -c 1.0 -s 0.05 -i examples/input/tree1_ -o examples/output/tree1_ -N 1 -S 240900 -P examples/input/tree1_distribHperP.txt

#Small example with host switch and no cospeciation
bin/cophylo.out -l 0.2 -m 0 -t 4.5 -c 0 -s 1.0 -i examples/input/example_ -o examples/output/host-switch_example_ -N 0 -S 3 -P examples/input/example_distribHperP.txt

#Bird-Malaria cophylogeny with random cospeciation probability and host switch rate
bin/cophylo.out -l 0.96 -m 0.59 -t 13.96 -i examples/input/bird-malaria_ -o examples/output/bird-malaria_ -N 50 -S 3 -P examples/input/bird-malaria_distribHperP.txt

#Run R script which plots the cophylogenies, creates and plots the Host-Parasite network, and computes the 51 network metrics used in the manuscript; requires packages ape and igraph
Rscript Rscript/compute.R  5 4.5 examples/input/example_ examples/output/cospeciation_example_ examples/cospeciation_example_
Rscript Rscript/compute.R  5 4.5 examples/input/example_ examples/output/host-switch_example_ examples/host-switch_example_
Rscript Rscript/compute.R  100 13.96 examples/input/bird-malaria_ examples/output/bird-malaria_ examples/bird-malaria_