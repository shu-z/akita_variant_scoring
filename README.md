Pipeline for generating mutated sequences for input into predictive models and its application using Akita to measure disruption to genome folding (Fudenberg et. al. 2020).
  
Install Akita and its dependencies:  
https://github.com/calico/basenji/tree/master/manuscripts/akita  
  
Packages needed:  
	- Numpy    
	- Pandas  
	- Scipy   
	- Tensorflow  
	- Basenji and its dependencies  
	- Cooltools  
	- Pysam   
	- Math  
	- Collection  
	- Bioseq  

<<<<<<< HEAD
Test run on a vcf of structural variants:  
> python score_var.py --in test/input/subset.somaticSV.vcf --shift_by -1 0 1 --file subset_SV --dir test/output



