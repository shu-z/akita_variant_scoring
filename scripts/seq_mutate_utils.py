#!/usr/bin/env python
# coding: utf-8



'''
Functions that accompany SuPreMo, to use after get_seq. 
Options to mutate the following after get_seq:
 - GC content of variant and flanking regions 
 - Shuffle nucleotide of variant and surrounding regions

'''

# # # # # # # # # # # # # # # # # # 
# # # # # Import packages # # # # #

import pandas as pd
import numpy as np
import random

import os
import io


# # # # # # # # # # # # # # # # # # 












# # # # # # # # # # # # # # # # # # 
# helper functions
# # # # # # # # # # # # # # # # # # 

def mutate_gc(seq, variant_start, variant_end, posflank, endflank, revcomp, mut_percent):
    """
    # note to self: this is different for BNDs -- can only change GC in flanking regions 
    # 
    """

    if not (0 <= mut_percent <= 100):
        raise ValueError("Mutation percentage must be between 0 and 100")
        
    if (variant_start-posflank)<0 or (variant_end+endflank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")
    
    mutated_sequence = list(seq)
    if posflank>0:
        mut_start=variant_start-posflank
        mut_end=variant_start
    if endflank>0:
        mut_start=variant_end
        mut_end=variant_end+endflank
        
    #right now, just assumes that revcomp does not deal with flanking sequences 
    #will need to fix it later
    if revcomp:
        #also make sure these are not one off 
        mut_start=len(seq)-variant_end
        mut_end=len(seq)-variant_start
        
    
    else:
        mut_start=variant_start
        mut_end=variant_end
    
    
    for i in range(mut_start, mut_end + 1):
        if random.randint(1, 100) <= mut_percent:
            current_nucleotide = mutated_sequence[i]
            if current_nucleotide in ['G', 'C']:
                mutated_sequence[i] = random.choice(['A', 'T'])
                
    return ''.join(mutated_sequence)

    
def gc_seq(seq, start, end):
    seq_substr=seq[start:(end+1)]
    gc_count = seq_substr.count('G') + seq_substr.count('C')
    return (gc_count / len(seq_substr)) 
    



def shuffle_nucs(seq, variant_start, variant_end, posflank, endflank, revcomp):
    if variant_start < 0 or variant_end >= len(seq) or variant_start > variant_end:
        raise ValueError("Invalid range of indices")
        
           
    if (variant_start-posflank)<0 or (variant_end+endflank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")
    
    mutated_sequence = list(seq)
    if posflank>0:
        shuff_start=variant_start-posflank
        shuff_end=variant_start
    if endflank>0:
        shuff_start=variant_end
        shuff_end=variant_end+endflank
        
    #right now, just assumes that revcomp does not deal with flanking sequences 
    #will need to fix it later
    if revcomp:
        #also make sure these are not one off 
        mut_start=len(seq)-variant_end
        mut_end=len(seq)-variant_start
    
    else:
        shuff_start=variant_start
        shuff_end=variant_end
          
        
    seq_list = list(seq)
    seq_to_shuffle = seq_list[shuff_start:shuff_end + 1]
    random.shuffle(seq_to_shuffle)
    seq_list[shuff_start:shuff_end + 1] = seq_to_shuffle

    return ''.join(seq_list)  