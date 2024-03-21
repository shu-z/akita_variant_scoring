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
import warnings


# # # # # # # # # # # # # # # # # # 


# this is a temporary fix 
warnings.filterwarnings("ignore")


def mutate_sequence(CHR, POS, END, SVTYPE, SVLEN, sequences_i, shift, revcomp, gc_mutate, shuffle_mutate, tf_mutate):
    """
    Applies sequence mutation function based on type of mutation specified in gc/shuffle/tf_mutate


    Arguments:
        gc_mutate: empty dictionary OR dictionary containing the parameters for mutating GC content near/around the variant
        shuffle_mutate: empty dictionary OR dictionary containing the parameters for shuffling nucleotides near/around the variant
        tf_mutate:empty dictionary OR dictionary containing the parameters for mutating TF binding motifs near/around the variant


    Returns:
        sequences_i: 


    """
    print('Mutating sequence');
    
    if gc_mutate:
        #we want to alter ref_seq in sequences_i
        #mutate_gc(sequences_i, up_flank, down_flank, revcomp, mut_percent):      
        #sequences_i_new=mutate_gc(sequences_i, gc_mutate['up_flank'], gc_mutate['down_flank'], revcomp, gc_mutate['mut_percent'])
        print('gc mutate', gc_mutate)
        sequences_i_new=mutate_gc(sequences_i, float(gc_mutate[0]), float(gc_mutate[1]), revcomp, int(gc_mutate[2]), SVLEN)

        
        #returns REF and GC-altered REF atm 
        return(sequences_i_new)
                        
                        

    if shuffle_mutate:
        #we want to alter ref_seq in sequences_i
        #shuffle_nucs(sequences_i, up_flank, down_flank, revcomp):      
        #sequences_i_new=shuffle_nucs(sequences_i, shuffle_mutate['up_flank'], shuffle_mutate['down_flank'], revcomp)
        sequences_i_new=shuffle_nucs(sequences_i, int(shuffle_mutate[0]), int(shuffle_mutate[1]), revcomp, SVLEN)

                        
        #returns REF and shuffled REF atm 
        return(sequences_i_new)

    
    if tf_mutate:
        #this is the argument 
        #tf_mutate={'up_flank': 1000, 'down_flank':1000, 'tf_mutate_type': 'INV', 'alter_variant': False, 'track_path':''}
        
        # sequences_i_new=mutate_TFs(CHR, POS, END, SVTYPE, SVLEN, sequences_i, revcomp, 
        #                            tf_mutate['tf_mutate_type'], tf_mutate['up_flank'], tf_mutate['down_flank'], 
        #                            tf_mutate['alter_variant'], tf_mutate['track_path'])

        sequences_i_new=mutate_TFs(CHR, POS, END, SVTYPE, SVLEN, sequences_i, revcomp, 
                                   tf_mutate[2], int(tf_mutate[0]), int(tf_mutate[1]), 
                                   tf_mutate[3], tf_mutate[4])
        return(sequences_i_new)




# # # # # # # # # # # # # # # # # # 
# helper functions
# # # # # # # # # # # # # # # # # # 

    
def gc_seq(seq, start, end):
    """
    This function returns the GC percentage of a sequence. Ignores bases with "N"
    """
    seq_substr=seq[start:(end+1)]
    gc_count = seq_substr.count('G') + seq_substr.count('C')
    n_count=seq_substr.count('N')
    return (gc_count / (len(seq_substr) - n_count) )



def mutate_gc(sequences_i, up_flank, down_flank, revcomp, mut_percent, SVLEN):
    """
    #note to self: this works on the REF seq, so the variant isn't altered at the moment 
    
    # note to self: this is different for BNDs -- can only change GC in flanking regions 
    # this does not current work for BNDs 


    Arguments:
        up_flank: if 
        down_flank:
        mut_percent: percentange of 
        

    Returns:
        seq_i
    """

    #note to self: need to remote these print statements later, they mess up the log 
    

    #pull out info from sequences_i
    seq=sequences_i[0]
    var_rel_pos=int(sequences_i[-1][0])
    variant_rel_end=var_rel_pos+float(SVLEN)    
    
    if not (0 <= mut_percent <= 100):
        raise ValueError("Mutation percentage must be between 0 and 100")
        
    if (var_rel_pos-up_flank)<0 or (variant_rel_end+down_flank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")


    
    mutated_sequence = list(seq)
    #print('mutated sequence pre', mutated_sequence[0:5])
    if up_flank>0:
        mut_start=var_rel_pos-up_flank
        mut_end=var_rel_pos
    if down_flank>0:
        mut_start=variant_rel_end
        mut_end=variant_rel_end+down_flank
        
    #right now, just assumes that revcomp does not deal with flanking sequences 
    #will need to fix it later
    if revcomp:
        #also make sure these are not one off 
        mut_start=len(seq)-variant_rel_end
        mut_end=len(seq)-var_rel_pos
        
    
    else:
        mut_start=var_rel_pos
        mut_end=variant_rel_end

    mut_start=int(mut_start)
    mut_end=int(mut_end)
    print('mut start end', mut_start, mut_end)
    print('old mut seq', mutated_sequence[mut_start:mut_end])
    
    print(f'Mutating gc at percentange {mut_percent}');
    for i in range(mut_start, mut_end + 1):
        if random.randint(1, 100) <= int(mut_percent):
            current_nucleotide = mutated_sequence[i]
            if current_nucleotide in ['G', 'C']:
                mutated_sequence[i] = random.choice(['A', 'T'])

    print('new mut seq', mutated_sequence[mut_start:mut_end])

    #join together new seq
    new_seq=''.join(mutated_sequence)
     
    #for now, replace the alt seq with the gc mutated seq
    sequences_i_list=list(sequences_i)
    sequences_i_list[1]=new_seq
    sequences_i_new=tuple(sequences_i_list)

                
    return sequences_i_new


    



def shuffle_nucs(sequences_i, up_flank, down_flank, revcomp, SVLEN):
    """
    This function shuffles the nucleotide order of a variant (and or the flanking regions up/downstream of the variant) 
    """

    #note to self: need to remote these print statements later, they mess up the log 
    print('Shuffling nucleotides...');


    #pull out info from sequences_i
    seq=sequences_i[0]
    var_rel_pos=sequences_i[-1][0]
    variant_rel_end=var_rel_pos+SVLEN


    
    if var_rel_pos < 0 or variant_rel_end >= len(seq) or var_rel_pos > variant_rel_end:
        raise ValueError("Invalid range of indices")
        
           
    if (var_rel_pos-up_flank)<0 or (variant_rel_end+down_flank)>len(seq):
        raise ValueError("Flanking indices exceed available sequence length")
    
    mutated_sequence = list(seq)
    if up_flank>0:
        shuff_start=var_rel_pos-up_flank
        shuff_end=var_rel_pos
    if down_flank>0:
        shuff_start=variant_rel_end
        shuff_end=variant_rel_end+down_flank
        
    #right now, just assumes that revcomp does not deal with flanking sequences 
    #will need to fix it later
    if revcomp:
        #also make sure these are not one off 
        shuff_start=len(seq)-variant_rel_end
        shuff_end=len(seq)-var_rel_pos
    
    else:
        shuff_start=var_rel_pos
        shuff_end=variant_rel_end

    shuff_start=int(shuff_start)
    shuff_end=int(shuff_end)
        
    seq_list = list(seq)
    seq_to_shuffle = seq_list[shuff_start:shuff_end + 1]
    random.shuffle(seq_to_shuffle)
    seq_list[shuff_start:shuff_end + 1] = seq_to_shuffle

    
    #join together new seq
    new_seq=''.join(seq_list)
                            
    #for now, replace the alt seq with the gc mutated seq
    sequences_i_list=list(sequences_i)
    sequences_i_list[1]=new_seq
    sequences_i_new=tuple(sequences_i_list)

                
    return sequences_i_new


################ motif alterations 

def shuffle_motifs(sequence, substr_ranges):
    """
    This function shuffles substrings of nucleotides within a genomic sequence

    
    """
    seq_list = list(sequence)

    # Shuffle nucleotides within each substring
    for start, end in substr_ranges:
        # Extract the substring
        substring = seq_list[start:end]

        # Shuffle the substring in-place
        random.shuffle(substring)

        # Update the sequence with the shuffled substring
        seq_list[start:end] = substring
    
    return ''.join(seq_list)



def invert_motifs(sequence, substr_ranges):
    """
    This function 
    """
    seq_list = list(sequence)  
    
    for start, end in substr_ranges:
        # Invert the substring within the list
        seq_list[start:end] = seq_list[start:end][::-1]
        
    # Convert the list back to a string
    return ''.join(seq_list)




def mutate_TFs(CHR, POS, END, SVTYPE, SVLEN, sequences_i, revcomp, tf_mutate_type, up_flank, down_flank, alter_variant, track_path):

    """
    This function mutates TFs within the variant and its flanking region
    

    Note to self: use a different function for shuffling motifs, that one is easier / doesn't need to count up to 1024s
    
    Arguments:
        up_flank: 
            - length of bp upstream of variant to look for TF motifs 
        down_flank:
            - length of bp downstream of variant to look for TF motifs    
        tf_mutate_type:
            - how we want to alter the tf motif
            -options: "INV", "shuffle"
        alter_variant: 
            - if True, we use the ALT sequence (variant included)
            - if False, we use the REF sequence (variant not included) 
            - THIS IS CURRENTLY NOT IMPLEMENETED
        track_path: 
            - path to TF binding track. For now, must be JASPAR-style bed tracks 

    Returns: 
        sequences_i: tuple containing:
            1. REF sequence (index 0)
            2. REF-altered sequence based on the TF motifs (index 1)
            3. variant POS indices relative to start for [REF, ALT] (index 2)

    
        Note to self: this is different for BNDs!!!!!!!!!!!!!!!!

        Note to self: might want to make this compatible with multiple TFs in the future?? maybe? would get complicated

    """

    if SVTYPE=='BND':
        #FIX THIS LATER 
        return(sequences_i)


    # check that if using flanking regions, they are within the length of the sequence    
    var_rel_pos=sequences_i[-1][0]     #this one is based on the REF sesquence 
    var_rel_end=var_rel_pos+SVLEN    #this one is based on the REF sesquence 
    sequence=sequences_i[0]

    if (var_rel_pos - up_flank) < 0 or (var_rel_end + down_flank) > len(sequence):
        raise ValueError("Flanking indices exceed available sequence length")

   
    # read in track 
    track_cols=['chr', 'start', 'end', 'tf', 'rel_score', '-log10(pval)', 'strand']
    track_df=pd.read_csv(track_path, sep='\t', header=None, skiprows=1, names=track_cols) 
    track_df['start']=track_df['start'].astype(int)
    track_df['end']=track_df['end'].astype(int)
    track_df['span']=track_df['end']-track_df['start']
    #tf_name=track_df.iloc[1,3]

    #get track df surrounding the variant
    print('track chr', track_df.chr[1])
    print("CHR", CHR)
    track_df_variant=track_df[(track_df.chr==CHR) & (track_df.start>(POS-up_flank)) & (track_df.end<(END+down_flank))]

    #later, may want to include option to only look at flanking regions and not include the variant itself 
    print('len track_df_variant', len(track_df_variant))
   
    #get locations of each motif for within each sequence  
    track_df_variant['diff_from_POS']=track_df_variant['start'] - POS
    track_df_variant['motif_start_idx']=var_rel_pos + track_df_variant['diff_from_POS']
    track_df_variant['motif_end_idx']=track_df_variant['motif_start_idx'] + track_df_variant['span']

    
    # if REVCOMP, need to invert these indices
    # DOUBLE CHECK THAT THIS IS CORRECT I AM WORRIED
    if revcomp:
        #also make sure these are not one off 
        track_df_variant['motif_start_idx_revcomp']=len(sequence)- track_df_variant['motif_end_idx']
        track_df_variant['motif_end_idx_revcomp']=len(sequence)- track_df_variant['motif_start_idx']
        
        #turn this into the list of tuples 
        tuple_coordinates = list(zip(track_df_variant['motif_start_idx_revcomp'], track_df_variant['motif_end_idx_revcomp']))

    else:
        #turn this into the list of tuples 
        tuple_coordinates = list(zip(track_df_variant['motif_start_idx'], track_df_variant['motif_end_idx']))



    #finally, alter sequences based on which alteration
    #current choices are only shuffling or inverting 


    if tf_mutate_type=='INV':
        new_seq=invert_motifs(sequence, tuple_coordinates)

    elif tf_mutate_type=='shuffle':
        new_seq=shuffle_motifs(sequence, tuple_coordinates)

    else:
        raise Exception("No suitable TF mutate option selected")

        
    #for now, replace the alt seq with the gc mutated seq
    sequences_i_list=list(sequences_i)
    sequences_i_list[1]=new_seq
    sequences_i=tuple(sequences_i_list)

    #returns REF and GC-altered REF atm 
    return(sequences_i)






    
    #the below only applies if we use deletions/duplications. for now, we don't need this 
    
    # #NEED TO CHECK THAT TOTAL MOTIF LENGTH IN REGION IS <1024
    # #IF NOT, NEED TO SAMPLE UNTIL MAX 1024 
    # if sum(track_df_variant['span']) > 100:

    #     average_span=np.mean(track_df_variant['span'])
    #     estimated_rows = min(int(100 / average_span), len(track_df_variant))

    #     # Sample rows from the DataFrame without replacement
    #     sampled_track_df = track_df_variant.sample(n=estimated_rows)

    #     # Calculate the total sum of spans in the sampled DataFrame
    #     total_span = sampled_track_df['span'].sum()

    #     # If the total span exceeds 100, remove rows until it's under 100
    #     while total_span > 100 and len(sampled_track_df) > 0:
            
    #         # get the row with the maximum span and remove 
    #         max_span_row = sampled_track_df.loc[sampled_track_df['span'].idxmax()]
    #         sampled_track_df = sampled_track_df.drop(max_span_row.name)
    
    #         # Update the total span
    #         total_span -= max_span_row['span']

    # else:
    #     sampled_track_df=track_df