Contains complex rearrangement data shared from Akdemir lab at MD Anderson, and code to analyze the rearrangement events with Akita.  


- data folder has complex SVs with the following files:
	- .fa files has raw fasta sequence ranging from <1Mb to a few Mb long 
	- .bed files has corresponding walk information for each .fa sequence
	- .mcool files for visualizing hi-c
	- Note that most of these complex connections involve a known cancer gene, which is indicated in the file name 

- bin folder has code for making akita predictions, visualizing hi-c, etc
	- variant calling pipeline adapted from Katie's code on github, https://github.com/shu-z/akita_variant_scoring
	- visualizations follow tutorials on cooltools, https://cooltools.readthedocs.io/en/latest/notebooks/viz.html

- results folder has ..... results .....


- figs folder has akita prediction plots  
