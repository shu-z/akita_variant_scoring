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



def mutate_sequence(CHR, POS, SVTYPE, SVLEN, sequences_i, shift, revcomp, gc_mutate, shuffle_mutate, tf_mutate) 
    """
    This function is calls other functions based on the type of sequence mutation 


    Arguments:
        gc_mutate: empty list or list containing
        shuffle_mutate: 
        tf_mutate:


    Returns:
        sequences_i: 


    """
    print('Mutating sequence');
    
    if gc_mutate:
        #we want to alter ref_seq in sequences_i
        #mutate_gc(seq, variant_start, variant_end, up_flank, down_flank, revcomp, mut_percent):
        var_rel_pos=sequences_i[-1][0]
        new_seq=seq_mutate_utils.mutate_gc(sequences_i[0], var_rel_pos, var_rel_pos+SVLEN, 
                                                           0,0, revcomp, 50)
                        
        #for now, replace the alt seq with the gc mutated seq
        sequences_i_list=list(sequences_i)
        sequences_i_list[1]=new_seq
        sequences_i=tuple(sequences_i_list)
                        
                        

    if shuffle_mutate:
        pass
                    
    if tf_mutate:
        var_rel_pos=sequences_i[-1][0]
        new_seq=seq_mutate_utils.shuffle_nucs(sequences_i[0], var_rel_pos, var_rel_pos+SVLEN, 
                                                              up_flank, down_flank, revcomp)
         
                        
        #for now, replace the alt seq with the gc mutated seq
        sequences_i[1]=new_seq 








# # # # # # # # # # # # # # # # # # 
# helper functions
# # # # # # # # # # # # # # # # # # 

    
def gc_seq(seq, start, end):
    """
    This function returns the GC percentage of a sequence. 
    """
    seq_substr=seq[start:(end+1)]
    gc_count = seq_substr.count('G') + seq_substr.count('C')
    return (gc_count / len(seq_substr)) 



def mutate_gc(seq, variant_start, variant_end, up_flank, down_flank, revcomp, mut_percent):
    """
    #note to self: this works on the REF seq, so the variant isn't altered at the moment 
    
    # note to self: this is different for BNDs -- can only change GC in flanking regions 
    # this does not current work for BNDs 

    #note to self: need to change the seq, variant_start, variant_end parameters to be taken from sequences_i

    Arguments:
        up_flank: if 
        down_flank:
        mut_percent: percentange of 
        

    Returns:
        seq_i
    """

    #note to self: need to remote these print statements later, they mess up the log 
    print('Mutating gc...');

    if not (0 <= mut_percent <= 100):
        raise ValueError("Mutation percentage must be between 0 and 100")
        
    if (variant_start-up_flank)<0 or (variant_end+down_flank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")
    
    mutated_sequence = list(seq)
    if up_flank>0:
        mut_start=variant_start-up_flank
        mut_end=variant_start
    if down_flank>0:
        mut_start=variant_end
        mut_end=variant_end+down_flank
        
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


    



def shuffle_nucs(seq, variant_start, variant_end, up_flank, down_flank, revcomp):
    """
    This function shuffles the nucleotide order of a variant (and or the flanking regions up/downstream of the variant) 
    """

    #note to self: need to remote these print statements later, they mess up the log 
    print('Shuffling nucleotides...');


    
    if variant_start < 0 or variant_end >= len(seq) or variant_start > variant_end:
        raise ValueError("Invalid range of indices")
        
           
    if (variant_start-up_flank)<0 or (variant_end+down_flank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")
    
    mutated_sequence = list(seq)
    if up_flank>0:
        shuff_start=variant_start-up_flank
        shuff_end=variant_start
    if down_flank>0:
        shuff_start=variant_end
        shuff_end=variant_end+down_flank
        
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