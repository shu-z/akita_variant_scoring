#! cd /pollard/home/shzhang/akita_variant_scoring/scripts/

import pandas as pd
import numpy as np
import os
from collections import Counter
import matplotlib.pyplot as plt
from scipy import stats
import sys
import math
from pybedtools import BedTool
import itertools
from matplotlib.backends.backend_pdf import PdfPages


import sys
os.chdir('/pollard/home/shzhang/akita_variant_scoring/')
sys.path.insert(0, '/pollard/home/shzhang/akita_variant_scoring/scripts/')

import plotting_utils
import get_Akita_scores_utils
import get_seq_utils

import matplotlib.gridspec as gridspec
from cooltools.lib.plotting import *
import itertools

import warnings
warnings.filterwarnings("ignore")
#warnings.filterwarnings("default")



seq_length=get_Akita_scores_utils.seq_length
bin_size=get_Akita_scores_utils.bin_size
target_length_cropped=get_Akita_scores_utils.target_length_cropped
half_patch_size=get_Akita_scores_utils.half_patch_size
map_length = bin_size * target_length_cropped


#for figures 
def save_figures_to_pdf(output_pdf_path, figures):
    with PdfPages(output_pdf_path) as pdf:
        for fig in figures:
            #plt.savefig("plot.pdf", bbox_inches='tight', pad_inches=0.5)
            pdf.savefig(fig, bbox_inches='tight', pad_inches=0.5)
            

            
#######################################################
#DATA AND SUCH 
            
# this is SVs from all CBTN tumors
# katie used 98.9% cutoff to define disruption, both mse and correlation
data_dir='/pollard/data/projects/CBTN_collab/'

CBTN_SV_scores=pd.read_csv(f'{data_dir}CBTN_SV_scores', sep='\t')
percentile_score = CBTN_SV_scores['mse_mean'].quantile(0.989)

chrom_lengths_path='/pollard/home/shzhang/akita_variant_scoring/data/chrom_lengths_hg38'
chrom_lengths=pd.read_table(chrom_lengths_path, header = None, names = ['CHROM', 'chrom_max'])
centromere_coords_path = f'/pollard/home/shzhang/akita_variant_scoring/data/centromere_coords_hg38'
centromere_coords = pd.read_table(centromere_coords_path, sep = '\t')


gene_annot_hg38 = pd.read_csv(f'/pollard/home/shzhang/akita_variant_scoring/data/gene_annot_hg38', sep = ',')
gene_annot_hg38 = gene_annot_hg38[['.' not in x for x in gene_annot_hg38.gene]]
gene_annot_hg38_BED = BedTool.from_dataframe(gene_annot_hg38[['chr', 'start', 'end', 'gene']])


#######################################################
#FUNCTIONS TO RUN ABC WEIGHTED SCORIGN


def get_variant_position(CHR, POS, var_len, half_left, half_right):
    
    '''
    Annotate variants based on position relative to chromosome arm ends (within half the size of the prediction window of end).
    One of 5 annotations: 
    1. chrom_mid: not near chromosome arm ends
    2. chrom_start: near chromosome start
    3. chrom_centro_left: near the left end of centromere (left arm)
    4. chrom_centro_right: near the right end of centromere (right arm)
    5. chrom_end: near chromosome end.
    5. centromere: overlaps centromere. These variants will not be processed
    
    '''
    
    # Define variant position with respect to chromosome start and end

    # Get last coordinate of chromosome
    chrom_max = int(chrom_lengths[chrom_lengths.CHROM == CHR[3:]]['chrom_max'].iloc[0]) 
    
    # Get centromere coordinate
    centro_start = int(centromere_coords[centromere_coords.chrom == CHR]['start'].iloc[0])
    centro_stop = int(centromere_coords[centromere_coords.chrom == CHR]['end'].iloc[0])

    # If variant too close to the beginning of chromosome
    if POS - half_left <= 0: 
        var_position = "chrom_start"
        
    # If variant in centromere
    elif POS >= centro_start and POS <= centro_stop:
        var_position = "centromere"
        
    # If variant too close to left end of centromere
    elif POS + var_len - 1 + half_right > centro_start and POS - half_left < centro_stop: 
        
        if POS > centro_stop:
            var_position = "chrom_centro_right"

        elif POS < centro_start:
            var_position = "chrom_centro_left"

    # If variant too close to the end of chromosome
    elif POS + var_len - 1 + half_right > chrom_max: 
        var_position = "chrom_end"
  
    else:
        var_position = "chrom_mid"
           
    return var_position


#altered from get_seq_utils.py
def get_BND_ALT_sense_coords(CHR_left, POS_left, left_position, CHR_right, POS_right, right_position, shift):
    
    '''
    This function applies to t[p[ and ]p]t BNDs.
    
    Get sequences corresponding to the left and right side of the alternate allele of the BND and two references.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''
    
    #print('running BND ALT sense')
    
    left_start = POS_left - half_patch_size + shift
    right_end = POS_right + half_patch_size + shift
    
    
    # Adjust start and end positions based on left and right variant positions
    adjust_BND_start_end_class = get_seq_utils.adjust_BND_start_end(CHR_left, POS_left, left_position, 
                                                      CHR_right, POS_right, right_position, 
                                                      left_start, right_end, shift)
    left_start, right_end = getattr(adjust_BND_start_end_class, 
                                    get_seq_utils.get_BND_ALT_sense_by_pos[left_position+'_'+right_position][0])()        
    
    
    # Get relative position of breakend with respect to start of map
    BND_rel_pos = POS_left - left_start
    
    if BND_rel_pos != seq_length/2:
        print("Warning: Variant not centered; too close to chromosome arm ends.")
    
    
    #ORIGINAL FUNCTION left and right sequences 
#     # Get sequences with adjusted sequence ends
#     ALT_left = fasta_open.fetch(CHR_left, left_start, POS_left).upper()
#     ALT_right = fasta_open.fetch(CHR_right, POS_right - 1, right_end - 1).upper() 


    BND_left_coord=(CHR_left, left_start, POS_left)
    BND_right_coord=(CHR_right, POS_right - 1, right_end - 1)
    
    return BND_left_coord, BND_right_coord, BND_rel_pos




#altered from get_seq_utils.py
def get_BND_ALT_antisense_left_coords(CHR_left, POS_left, left_position, CHR_right, POS_right, right_position, shift):      
    
    '''
    This function applies to [p[t BNDs.
    
    Get sequences corresponding to the left and right side of the alternate allele of the BND and two references.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''
    
    #print('running BND ALT antisense left coords')


    left_start = POS_left - (- half_patch_size + shift)
    right_end = POS_right + half_patch_size + shift

    
    # Adjust start and end positions based on left and right variant positions
    adjust_BND_start_end_class = get_seq_utils.adjust_BND_start_end(CHR_left, POS_left, left_position, 
                                                      CHR_right, POS_right, right_position, 
                                                      left_start, right_end, shift)
    left_start, right_end = getattr(adjust_BND_start_end_class, 
                                    get_seq_utils.get_BND_ALT_sense_by_pos[left_position+'_'+right_position][1])()
    
    
    # Get relative position of breakend with respect to start of map
    BND_rel_pos = left_start - POS_left
    
    if BND_rel_pos != seq_length/2:
        print("Warning: Variant not centered; too close to chromosome arm ends.")
      
        
#     # Get sequences with adjusted sequence ends
#     ALT_left_revcomp = fasta_open.fetch(CHR_left, POS_left - 1, left_start - 1).upper()
#     ALT_left = str(Seq(ALT_left_revcomp).reverse_complement())

#     ALT_right = fasta_open.fetch(CHR_right, POS_right - 1, right_end - 1).upper()


    #left coord is antisense!
    BND_left_coord=(CHR_left, POS_left - 1,  left_start - 1)
    BND_right_coord=(CHR_right, POS_right - 1, right_end - 1)
    
    return BND_left_coord, BND_right_coord, BND_rel_pos


#altered from get_seq_utils.py
def get_BND_ALT_antisense_right_coords(CHR_left, POS_left, left_position, CHR_right, POS_right, right_position, shift): 
    
    '''
    This function applies to t]p] BNDs.
    
    Get sequences corresponding to the left and right side of the alternate allele of the BND and two references.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''
    
    #print('running BND ALT antisense right')


    left_start = POS_left - half_patch_size + shift
    right_end = POS_right - (half_patch_size + shift)

 
    # Adjust start and end positions based on left and right variant positions
    adjust_BND_start_end_class = get_seq_utils.adjust_BND_start_end(CHR_left, POS_left, left_position, 
                                                      CHR_right, POS_right, right_position, 
                                                      left_start, right_end, shift)
    left_start, right_end = getattr(adjust_BND_start_end_class, 
                                    get_seq_utils.get_BND_ALT_sense_by_pos[left_position+'_'+right_position][2])()
    
    
    # Get relative position of breakend with respect to start of map
    BND_rel_pos = POS_left - left_start
    
    if BND_rel_pos != seq_length/2:
        print("Warning: Variant not centered; too close to chromosome arm ends.")
      

    # Get sequences with adjusted sequence ends
#     ALT_left = fasta_open.fetch(CHR_left, left_start, POS_left).upper()
#     ALT_right_revcomp = fasta_open.fetch(CHR_right, right_end, POS_right).upper()
#     ALT_right = str(Seq(ALT_right_revcomp).reverse_complement())

            
    #right coord is antisense!
    BND_left_coord=(CHR_left, left_start, POS_left)
    BND_right_coord=(CHR_right, right_end, POS_right)

    
    return BND_left_coord, BND_right_coord, BND_rel_pos
            


    
    
    
#this is taken from get_seq_utils.py
#ask Katie G -- in order to account for BND roi positioning, need to get start coordinates for both sides of BND
#currently only saves rel_pos_map for the left side 

def get_coords_BND(CHR, POS, REF, ALT, shift, revcomp):
    
    '''
    shu edit 
    MODFIED FROM get_sequences_BND
    USED with weighting 
    
    This function applies all BNDs.
    
    Get sequences corresponding to the alternate BND allele (joined breakends) and two references (one from each locus that was joined). To do so, get info for the second region (CHR2, POS2), get sequences, input the alternate allele at the breakend, and crop the ends of the sequence if necessary.
    Also get BND_rel_pos which is the relative position of the breakend in the generated reference sequences. 
    That corresponds to the end of the left alternate sequnce and the beginning of the right alternate sequence. 
    
    Note: 1 is subtracted in fasta_open.fetch because the hg38 fasta file is 0-based and all the coordinates are 1-based.
    If it's not subtracted, it's so POS is included in the sequence if it's at the end

    '''
    
    var_position = get_variant_position(CHR, POS, 0, half_patch_size - shift, half_patch_size + shift)
    
    if var_position == 'centromere':
        raise ValueError('Centromeric variant.')

    if '[' in ALT:
        
        CHR2 = ALT.split(':')[0].split('[')[1]
        POS2 = int(ALT.split('[')[1].split(':')[1])
        
        var_position2 = get_variant_position(CHR2, POS2, 0, half_patch_size - shift, half_patch_size + shift)
        
        if var_position2 == 'centromere':
            raise ValueError('Centromeric variant.')

        if ALT[0] == '[':
            print('BND is [p[t, antisense is left')
            #  [p[t
            antisense='left'
            ALT_t = ALT.split('[')[2]
            BND_left_coord, BND_right_coord, BND_rel_pos = get_BND_ALT_antisense_left_coords(CHR2, POS2, var_position2, 
                                                                                                      CHR, POS, var_position, 
                                                                                                      shift)

        elif ALT[0] != '[':
            print('BND is t[p[, antisense is none')
            # t[p[
            antisense='none'
            ALT_t = ALT.split('[')[0]
            BND_left_coord, BND_right_coord, BND_rel_pos = get_BND_ALT_sense_coords(CHR, POS, var_position, 
                                                                                             CHR2, POS2, var_position2, 
                                                                                             shift)


    elif ']' in ALT:
        
        CHR2 = ALT.split(':')[0].split(']')[1]
        POS2 = int(ALT.split(']')[1].split(':')[1])
        
        var_position2 = get_variant_position(CHR2, POS2, 0, half_patch_size - shift, half_patch_size + shift)
        
        if var_position2 == 'centromere':
            raise ValueError('Centromeric variant.')

        if ALT[0] == ']':
            print('BND is ]p]t, antisense is none')

            # ]p]t
            antisense='none'
            ALT_t = ALT.split(']')[2]
            BND_left_coord, BND_right_coord, BND_rel_pos = get_BND_ALT_sense_coords(CHR2, POS2, var_position2, 
                                                                                             CHR, POS, var_position, 
                                                                                             shift)
            
        elif ALT[0] != ']':
            print('BND is t]p], antisense is right')

            # t]p]            
            antisense='right'
            ALT_t = ALT.split(']')[0]
            BND_left_coord, BND_right_coord, BND_rel_pos = get_BND_ALT_antisense_right_coords(CHR, POS, var_position, 
                                                                                              CHR2, POS2, var_position2, 
                                                                                                  shift)
    #print('ALT_t', ALT_t)
  
#     if ALT_t != REF:
        
#         # Since we are already including REF, remove it from ALT_t
#         if ALT_t.startswith(REF):
#             ALT_t = ALT_t[1:]
#         elif ALT_t.endswith(REF):
#             ALT_t = ALT_t[:-1]
#         else:
#             raise ValueError('Unexpected format: BND ALT does not include REF.')
    
#         # Create ALT sequence
#         ALT_seq = ALT_left + ALT_t + ALT_right
        
#         # chop off the sides if longer than seq_length
#         ALT_seq = crop_sequence(ALT_seq, seq_length)
        
#     else:
#         # Create ALT sequence
#         ALT_seq = ALT_left + ALT_right
        
        
#     if revcomp:
#         REF_for_left, REF_for_right, ALT_seq = [str(Seq(x).reverse_complement()) for x in [REF_for_left, REF_for_right, ALT_seq]]

        
    #check length of BND left and right coords
    #print('bnd left coord length', BND_left_coord[2]-BND_left_coord[1])
    #print('bnd right coord length', BND_right_coord[2]-BND_right_coord[1])
    #print('bnd sum length')
                                                                                              
                                                                                              
    return BND_left_coord, BND_right_coord, BND_rel_pos, antisense



#edited input to include roi_coords_BED
#10/14/24 editing to also account for BNDs 
def get_roi_in_map_BND(coord_args, rel_pos_map, roi_coords_BED):
    
    
    '''
    Get bins for regions of interest (roi) that fall in the predicted map for BND Variants.
    
    Requires getting coordinates for two sides/coordinates of BND
    
    '''
    
#     this is not necessary, map start coord gives coordinates of only left side, rel_pos_map is same as BND_rel_pos
#     map_start_coord = map_start_coord[0]
#     rel_pos_map = rel_pos_map[1]
    
        
    #get BND coordinates 
    BND_left_coord, BND_right_coord, BND_rel_pos, antisense = get_coords_BND(*coord_args)
    #print(BND_rel_pos)

    #combine two map regions into dataframe
    map_regions_df=pd.DataFrame({'CHR' : [BND_left_coord[0], BND_right_coord[0]],
                                 'start' : [BND_left_coord[1], BND_right_coord[1]],
                                 'end' : [BND_left_coord[2], BND_right_coord[2]]})

    #turn into BED 
    map_region_BED = BedTool.from_dataframe(map_regions_df)

        

    #intersect roi coords with combined BND BED
    roi_in_map_BED = map_region_BED.intersect(roi_coords_BED, wa = True, wb = True)

    #nothing in map 
    if roi_in_map_BED == '':
        roi_in_map = ''

    else:
        #edit to include abc score 
        roi_in_map = (roi_in_map_BED
                        .to_dataframe()
                      #adding chrom 
                        .rename(columns = {'name': 'Chrom', 'score': 'Start', 'strand':'End', 
                                           'thickStart':'roi_id', 'thickEnd': "ABC.score"})[['Chrom', 'Start', 'End', 'roi_id', 'ABC.score']])
                     
                    
        
        #split into left and right sides 
        roi_in_map_left=roi_in_map[roi_in_map.Chrom==BND_left_coord[0]]
        roi_in_map_right=roi_in_map[roi_in_map.Chrom==BND_right_coord[0]]
        
        #get bin coordinates with respect to map for ROI regions 
        #THIS NEEDS TO BE FLIPPED BASED ON SENSE OR ANTISENSE ON LEFT/RIGHT
        #do for left and right separately                 
        if len(roi_in_map_left)>0:
            
            if antisense=='left':
                #flip
                
                #'BND is [p[t

                #get relative bin start position, with respect to the END of BND_left_coord
                
                #need to make this otherwise start gets overwritten, and need to call it again 
                roi_left_orig_start=roi_in_map_left.Start
                roi_left_orig_end=roi_in_map_left.End
                
                roi_in_map_left.Start = [math.ceil(x/bin_size) for x in (BND_left_coord[2] - roi_left_orig_end)]
                roi_in_map_left.loc[roi_in_map_left.Start < 0,'Start'] = 0

                
                #get relative bin end position, which is really the roi start 
                roi_in_map_left.End = [math.ceil(x/bin_size) for x in (BND_left_coord[2] - roi_left_orig_start)]
                
                
                #remove starts if after end of left half of chr 
                roi_in_map_left = roi_in_map_left[roi_in_map_left.Start <= rel_pos_map-1]
                
                #remove if it takes up past where the rel pos is 
                roi_in_map_left.loc[roi_in_map_left.End > rel_pos_map-1,'End'] = rel_pos_map-1
                
            else:
                #get relative bin start position
                roi_in_map_left.Start = [math.ceil(x/bin_size) for x in (roi_in_map_left.Start - BND_left_coord[1])]
                roi_in_map_left.loc[roi_in_map_left.Start < 0,'Start'] = 0
                #remove starts if after end of left half of chr 
                roi_in_map_left = roi_in_map_left[roi_in_map_left.Start <= rel_pos_map-1]
                
                #get relative bin end position
                roi_in_map_left.End = [math.ceil(x/bin_size) for x in (roi_in_map_left.End - BND_left_coord[1])]
                #remove if it takes up past where the rel pos is 
                roi_in_map_left.loc[roi_in_map_left.End > rel_pos_map-1,'End'] = rel_pos_map-1

        else:
            roi_in_map_left = pd.DataFrame(columns=['Chrom', 'Start', 'End', 'roi_id', 'ABC.score'])
        
        #print('roi in map left,', roi_in_map_left)

        
        #check right side map         
        if len(roi_in_map_right)>0:
            if antisense=='right':
                #flip
                #BND is t]p]
                
                roi_right_orig_start=roi_in_map_right.Start
                roi_right_orig_end=roi_in_map_right.End
                
                

                #get relative bin start position, with respect to END of BND_right_coord 
                roi_in_map_right.Start = [math.ceil(x/bin_size) for x in (BND_right_coord[2]-roi_right_orig_end)]
                #add based on where right starts (if middle, is 224)
                roi_in_map_right.Start = roi_in_map_right.Start + rel_pos_map
                #if start extends before right half of map, make it start right at BND_rel_pos
                roi_in_map_right.loc[roi_in_map_right.Start < rel_pos_map,'Start'] = rel_pos_map

                
                #get relative bin end position
                roi_in_map_right.End = [math.ceil(x/bin_size) for x in (BND_right_coord[2]-roi_right_orig_start)]
                roi_in_map_right.End = roi_in_map_right.End + rel_pos_map

                #crop end if it takes up past where the rel pos is 
                roi_in_map_right.loc[roi_in_map_right.End > target_length_cropped-1,'End'] = target_length_cropped-1

                
                #remove starts if after end of map 
                roi_in_map_right = roi_in_map_right[roi_in_map_right.Start <=target_length_cropped-1]
                
                
            else: 

                #get relative bin start position 
                roi_in_map_right.Start = [math.ceil(x/bin_size) for x in (roi_in_map_right.Start - BND_right_coord[1])]
                #add based on where right starts (if middle, is 224)
                roi_in_map_right.Start = roi_in_map_right.Start + rel_pos_map
                #if start extends before right half of map, make it start right at BND_rel_pos
                roi_in_map_right.loc[roi_in_map_right.Start < rel_pos_map,'Start'] = rel_pos_map
                #remove starts if after end of map 
                roi_in_map_right = roi_in_map_right[roi_in_map_right.Start <=target_length_cropped-1]
                
                #get relative bin end position
                roi_in_map_right.End = [math.ceil(x/bin_size) for x in (roi_in_map_right.End - BND_right_coord[1])]
                roi_in_map_right.End = roi_in_map_right.End + rel_pos_map

                #crop end if it takes up past where the rel pos is 
                roi_in_map_right.loc[roi_in_map_right.End > target_length_cropped-1,'End'] = target_length_cropped-1

        
        
        else:
            roi_in_map_right = pd.DataFrame(columns=['Chrom', 'Start', 'End','roi_id', 'ABC.score'])
            
            
        #print('roi in map right', roi_in_map_right)

                                                                     
                                                                    
    #combine roi in maps again 
    roi_in_map=pd.concat([roi_in_map_left, roi_in_map_right], axis=0)
    #if empty dataframe one side will mess things up s
    roi_in_map['Start'] = roi_in_map['Start'].astype(int)
    roi_in_map['End'] = roi_in_map['End'].astype(int)

                                                                       
    #calculate widths 
    roi_in_map['width'] = [x-y+1 for x,y in zip(roi_in_map.End, roi_in_map.Start)]
    roi_in_map.reset_index(drop=True, inplace = True)


    #print('roi in map', roi_in_map)
    return roi_in_map




def get_genes_in_map_BND(coord_args, rel_pos_map, genome = 'hg38'):
    
    
    '''
    Get protein coding genes in map region.
    
    '''

    #get BND coordinates 
    BND_left_coord, BND_right_coord, BND_rel_pos, antisense = get_coords_BND(*coord_args)

    #combine two map regions into dataframe
    map_regions_df=pd.DataFrame({'CHR' : [BND_left_coord[0], BND_right_coord[0]],
                                 'start' : [BND_left_coord[1], BND_right_coord[1]],
                                 'end' : [BND_left_coord[2], BND_right_coord[2]]})

    #turn into BED 
    map_region_BED = BedTool.from_dataframe(map_regions_df)
    
    
    if genome == 'hg38':
        gene_annot_BED = gene_annot_hg38_BED
    elif genome == 'hg19':
        gene_annot_BED = gene_annot_hg19_BED
        
    genes_in_map_BED = map_region_BED.intersect(gene_annot_BED, wa = True, wb = True)

    
    if genes_in_map_BED == '':
        genes_in_map = pd.DataFrame(columns=['Start', 'End', 'Gene', 'width'])
        
    else:
        genes_in_map = (genes_in_map_BED.to_dataframe()
                        .rename(columns = {'name': 'Chrom', 'score': 'Start', 'strand':'End', 'thickStart':'Gene'})
                        [['Chrom', 'Start', 'End', 'Gene']])
                
        #split into left and right sides 
        genes_in_map_left=genes_in_map[genes_in_map.Chrom==BND_left_coord[0]]
        genes_in_map_right=genes_in_map[genes_in_map.Chrom==BND_right_coord[0]]
        
        #get bin coordinates with respect to map for ROI regions 
        #THIS NEEDS TO BE FLIPPED BASED ON SENSE OR ANTISENSE ON LEFT/RIGHT
        #do for left and right separately                 
        if len(genes_in_map_left)>0:
            
            if antisense=='left':
                #flip
                
                #'BND is [p[t

                #get relative bin start position, with respect to the END of BND_left_coord
                
                #need to make this otherwise start gets overwritten, and need to call it again 
                gene_left_orig_start=genes_in_map_left.Start
                gene_left_orig_end=genes_in_map_left.End
                
                genes_in_map_left.Start = [math.ceil(x/bin_size) for x in (BND_left_coord[2] - gene_left_orig_end)]
                genes_in_map_left.loc[genes_in_map_left.Start < 0,'Start'] = 0

                
                #get relative bin end position, which is really the roi start 
                genes_in_map_left.End = [math.ceil(x/bin_size) for x in (BND_left_coord[2] - gene_left_orig_start)]
                

                
                
                #remove starts if after end of left half of chr 
                genes_in_map_left = genes_in_map_left[genes_in_map_left.Start <= rel_pos_map-1]
                
                #remove if it takes up past where the rel pos is 
                genes_in_map_left.loc[genes_in_map_left.End > rel_pos_map-1,'End'] = rel_pos_map-1
                
                genes_in_map_left['antisense']='antisense'

                
            else:
                #get relative bin start position
                genes_in_map_left.Start = [math.ceil(x/bin_size) for x in (genes_in_map_left.Start - BND_left_coord[1])]
                genes_in_map_left.loc[genes_in_map_left.Start < 0,'Start'] = 0
                #remove starts if after end of left half of chr 
                genes_in_map_left = genes_in_map_left[genes_in_map_left.Start <= rel_pos_map-1]
                
                #get relative bin end position
                genes_in_map_left.End = [math.ceil(x/bin_size) for x in (genes_in_map_left.End - BND_left_coord[1])]
                #remove if it takes up past where the rel pos is 
                genes_in_map_left.loc[genes_in_map_left.End > rel_pos_map-1,'End'] = rel_pos_map-1
                
                genes_in_map_left['antisense']='sense'


        else:
            genes_in_map_left = pd.DataFrame(columns=['Chrom', 'Start', 'End', 'Gene'])
        
        
        #check right side map         
        if len(genes_in_map_right)>0:
            if antisense=='right':
                #flip
                #BND is t]p]
                
                #print('BND right coord', BND_right_coord)
                
                #print('genes in map right', genes_in_map_right)
                
                genes_right_orig_start=genes_in_map_right.Start
                genes_right_orig_end=genes_in_map_right.End

                #get relative bin start position, with respect to END of BND_right_coord 
                genes_in_map_right.Start = [math.ceil(x/bin_size) for x in (BND_right_coord[2]-genes_right_orig_end)]
                #add based on where right starts (if middle, is 224)
                genes_in_map_right.Start = genes_in_map_right.Start + rel_pos_map
                #if start extends before right half of map, make it start right at BND_rel_pos
                genes_in_map_right.loc[genes_in_map_right.Start < rel_pos_map,'Start'] = rel_pos_map

                #print('genes in map right', genes_in_map_right)

                
                #get relative bin end position
                genes_in_map_right.End = [math.ceil(x/bin_size) for x in (BND_right_coord[2]-genes_right_orig_start)]
                genes_in_map_right.End = genes_in_map_right.End + rel_pos_map

                #crop end if it takes up past where the rel pos is 
                genes_in_map_right.loc[genes_in_map_right.End > target_length_cropped-1,'End'] = target_length_cropped-1

                
                #remove starts if after end of map 
                genes_in_map_right = genes_in_map_right[genes_in_map_right.Start <=target_length_cropped-1]
                
                #print('genes in map right', genes_in_map_right)

                
                genes_in_map_right['antisense']='antisense'
                
                
            else: 

                #get relative bin start position 
                genes_in_map_right.Start = [math.ceil(x/bin_size) for x in (genes_in_map_right.Start - BND_right_coord[1])]
                #add based on where right starts (if middle, is 224)
                genes_in_map_right.Start = genes_in_map_right.Start + rel_pos_map
                #if start extends before right half of map, make it start right at BND_rel_pos
                genes_in_map_right.loc[genes_in_map_right.Start < rel_pos_map,'Start'] = rel_pos_map
                #remove starts if after end of map 
                genes_in_map_right = genes_in_map_right[genes_in_map_right.Start <=target_length_cropped-1]
                
                #get relative bin end position
                genes_in_map_right.End = [math.ceil(x/bin_size) for x in (genes_in_map_right.End - BND_right_coord[1])]
                genes_in_map_right.End = genes_in_map_right.End + rel_pos_map

                #crop end if it takes up past where the rel pos is 
                genes_in_map_right.loc[genes_in_map_right.End > target_length_cropped-1,'End'] = target_length_cropped-1

                genes_in_map_right['antisense']='sense'

        
        else:
            genes_in_map_right = pd.DataFrame(columns=['Chrom', 'Start', 'End','Gene', 'antisense'])
            
            

                                                                     
                                                                    
    #combine roi in maps again 
    genes_in_map=pd.concat([genes_in_map_left, genes_in_map_right], axis=0)
    #if empty dataframe one side will mess things up s
    genes_in_map['Start'] = genes_in_map['Start'].astype(int)
    genes_in_map['End'] = genes_in_map['End'].astype(int)

                                                                       
    #calculate widths 
    genes_in_map['width'] = [x-y+1 for x,y in zip(genes_in_map.End, genes_in_map.Start)]
    genes_in_map.reset_index(drop=True, inplace = True)

            
            
    #print('gnees in map', genes_in_map)
    return genes_in_map




def get_weighted_ABC_score_BND(disruption_track, roi_in_map):

    '''
    I THINK THIS STAYS THE SAME, REGARDLESS OF ABC SCORE WEIGHTING OR NOT BUT NEED TO CHECK BACK IN ON THIS 
    '''

    #print('roi in map', roi_in_map)
    # Get the bins that correspond regions of interest (roi)
    
    roi_bins = []
    abc_bins= []
    for i in range(len(roi_in_map)):
        
        bins_i = list(range(roi_in_map.iloc[i].Start, roi_in_map.iloc[i].End + 1))

        abc_i= [roi_in_map.iloc[i]['ABC.score']] * len(bins_i)
        roi_bins.append(bins_i)
        abc_bins.append(abc_i)
        

    #turn into df
    roi_binlist=[x for sublist in roi_bins for x in sublist]
    abc_binlist=[x for sublist in abc_bins for x in sublist]
    bin_score_df=pd.DataFrame({'bin_num': roi_binlist, 'abc_score': abc_binlist})
    
    
    #get average ABC Score per bin, then subset to unique bins 
    bin_score_df['avg_abc'] = bin_score_df.groupby('bin_num')['abc_score'].transform('mean')
    bin_score_df_unique=bin_score_df.drop_duplicates(subset=('bin_num', 'avg_abc'))
    
    #initiate weight track of zeros (no enhancer, not included in score )
    background_weight = 0
    weight_track = np.array([background_weight]*len(disruption_track), dtype=float)
    
    #loop through bins, and abc score to weights track
    for row in bin_score_df_unique.itertuples():
        weight_track[row.bin_num] = row.avg_abc
    #print('weight track', weight_track[25:30])
        
    # Add nans to weight track
    weight_track[np.isnan(disruption_track)] = 0

    if len(roi_in_map) == 0:
        roi_id_values = ['']
    else:
        roi_id_values = roi_in_map.roi_id.values

    # Mask nan values in disruption track
    disruption_track_masked = np.ma.MaskedArray(disruption_track, mask=np.isnan(disruption_track))

    return np.ma.average(disruption_track_masked, weights=weight_track), weight_track, roi_id_values




def plot_maps_abc_BND(maps, genes_in_map, lines, roi_in_map, disruption_track, weight_track, title, scale = 1, gene_rows = 3):
    
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. 
    For duplications, there will be 3 lines marking the two regions that are duplicates.    
    Plot genes that match the regions in the reference map.
    Plot disruption score tracks.
    
    '''

    if len(genes_in_map)>0:
        gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    
    if len(roi_in_map)>0:
        roi_track = list(roi_in_map[['Start', 'width']].to_records(index = False))
    
    n = maps[0].shape[0]

    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    #plot_width1D determines height of the disruption, activity, etc tracks 
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width, .15, 
                                               plot_width1D, .15, plot_width1D, .15, plot_width1D, .15, plot_width1D/2])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        #map_plot=mat_from_vector(maps[i])
        #plot_upright(plt, map_plot, lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        plotting_utils.pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)

     
    #gets difference map 
    #print('plotting difference ')
    plt.subplot(gs[4,0])
    plotting_utils.pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
    #map_diff_plot=mat_from_vector(maps[0] - maps[1])

    #plot_upright(plt, map_diff_plot, lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)



    # Disruption track
   # print ('plotting disruption')
    plt.subplot(gs[6,0])
    plt.plot(list(range(target_length_cropped)), disruption_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.ylabel('Disruption') 

    # Weight track
    #print('plotting weights')
    plt.subplot(gs[8,0])
    plt.plot(list(range(target_length_cropped)), weight_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.ylabel('ABC Activity') 

    #plt.yticks([])
    

    # Gene track
    if len(genes_in_map)>0:
        ax1 = fig.add_subplot(gs[10,0])


        gene_bar_height = 15*gene_rows/3
        bar_locations = []

        for bar_level in range(gene_rows):

            bar_location = (bar_level+1)*(gene_bar_height/gene_rows) - 0.5
            ax1.broken_barh(gene_track[bar_level:][::gene_rows], 
                            (bar_location, 0.5*scale), 
                            facecolors='tab:blue')
            bar_locations.append(bar_location)


        bar_locations = (bar_locations*math.ceil(len(genes_in_map)/gene_rows))[:len(genes_in_map)]

        for i in range(len(genes_in_map)):
            bar_location = bar_locations[i]
            gene = genes_in_map.loc[i,'Gene']
            location = genes_in_map.loc[i,'Start']
            
            annotation_color = 'red' if genes_in_map.loc[i, 'antisense'] == 'antisense' else 'black'


            ax1.annotate(gene, (location,bar_location), # annotate gene at the start
                         rotation = 45, 
                         ha = "right", va = "top", # set horizontal and vertical alignment
                         annotation_clip = False, # keep annotation if outside of window
                         fontsize = 8*scale, color=annotation_color) 
            
        if lines is not None:
            for line in lines:
                plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
        plt.ylabel('Genes',rotation=90)
        plt.ylim([0,gene_bar_height])
        plt.xlim([0,448])
        plt.axis('off')
        
    
    fig.suptitle(title, fontsize=12, y=1.15)


    #plt.show()
    return(fig)




#edited input to include roi_coords_BED
#10/14/24 editing to also account for BNDs 
def get_roi_in_map(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN, roi_coords_BED):
    
    
    '''
    Get bins for regions of interest (roi) that fall in the predicted map.
    
    '''

    #shu commenting out because this doesn't seem necessary? but could be wrong
    #jk this is because of changing from int to list situation 
    #For duplications, we want to use regions that correspond to the masked REF matrix
    if SVTYPE == 'DUP':
        map_start_coord = map_start_coord[1]
        rel_pos_map = rel_pos_map[0]
    else:
        map_start_coord = map_start_coord[0]
        rel_pos_map = rel_pos_map[1]
        
    
    
    
    map_region_BED = BedTool.from_dataframe(pd.DataFrame({'CHR' : [CHR],
                                                          'start' : [map_start_coord],
                                                          'end' : [map_start_coord + seq_length]}))

    roi_in_map_BED = map_region_BED.intersect(roi_coords_BED, wa = True, wb = True)

    if roi_in_map_BED == '':
        roi_in_map = ''

    else:
        #edit to include abc score 
        roi_in_map = (roi_in_map_BED
                        .to_dataframe()
                        .rename(columns = {'score': 'Start', 'strand':'End', 
                                           'thickStart':'roi_id', 'thickEnd': "ABC.score"})[['Start', 'End', 'roi_id', 'ABC.score']])
                     
                      
        #print(roi_in_map)
        
        roi_in_map.Start = [math.ceil(x/bin_size) for x in (roi_in_map.Start - map_start_coord)]
        roi_in_map.loc[roi_in_map.Start < 0,'Start'] = 0
        roi_in_map = roi_in_map[roi_in_map.Start <= target_length_cropped-1]
        roi_in_map.End = [math.ceil(x/bin_size) for x in (roi_in_map.End - map_start_coord)]
        
        
        if SVTYPE == 'DUP' and abs(int(SVLEN)) > bin_size/2:

            roi_in_map.loc[(roi_in_map.Start >= rel_pos_map) & 
                              (roi_in_map.End > rel_pos_map),
                              ['Start', 'End']] = roi_in_map.loc[(roi_in_map.Start >= rel_pos_map) & 
                                                                  (roi_in_map.End > rel_pos_map),
                                                                  ['Start', 'End']] + math.ceil(SVLEN/bin_size)


            middle_rois = roi_in_map.loc[(roi_in_map.Start < rel_pos_map) & 
                              (roi_in_map.End > rel_pos_map)].roi_id


            for middle_roi in middle_rois:

                left = roi_in_map.loc[roi_in_map.roi_id == middle_roi]
                left.End = rel_pos_map
                right = roi_in_map.loc[roi_in_map.roi_id == middle_roi]
                right.Start = rel_pos_map + math.ceil(SVLEN/bin_size)
                right.End = right.End + math.ceil(SVLEN/bin_size)

                roi_in_map = pd.concat([roi_in_map[roi_in_map.roi_id != middle_roi], left, right], axis = 0).reset_index(drop=True)
            
            # Remove genes that left prediction window
            roi_in_map = roi_in_map[~((roi_in_map.End > target_length_cropped) &
                                        (roi_in_map.Start > target_length_cropped))]
            
            # Crop end of genes that left prediction window
            roi_in_map.loc[(roi_in_map.End > target_length_cropped),'End'] = target_length_cropped

            roi_in_map = roi_in_map[~np.isnan(roi_in_map.Start)]
            roi_in_map.Start = [int(x) for x in roi_in_map.Start]
            roi_in_map = roi_in_map[~np.isnan(roi_in_map.End)]
            roi_in_map.End = [int(x) for x in roi_in_map.End]
            
            roi_in_map.reset_index(drop=True, inplace = True)
        
        roi_in_map.loc[roi_in_map.End > target_length_cropped-1,'End'] = target_length_cropped-1
        roi_in_map['width'] = [x-y+1 for x,y in zip(roi_in_map.End, roi_in_map.Start)]
            
            
    return roi_in_map




def get_weighted_ABC_score(disruption_track, roi_in_map):


    #print('roi in map', roi_in_map)
    # Get the bins that correspond regions of interest (roi)
    
    roi_bins = []
    abc_bins= []
    for i in range(len(roi_in_map)):
        
        bins_i = list(range(roi_in_map.iloc[i].Start, roi_in_map.iloc[i].End + 1))

        abc_i= [roi_in_map.iloc[i]['ABC.score']] * len(bins_i)
        roi_bins.append(bins_i)
        abc_bins.append(abc_i)
        

    #turn into df
    roi_binlist=[x for sublist in roi_bins for x in sublist]
    abc_binlist=[x for sublist in abc_bins for x in sublist]
    bin_score_df=pd.DataFrame({'bin_num': roi_binlist, 'abc_score': abc_binlist})
    
    
    #get average ABC Score per bin, then subset to unique bins 
    bin_score_df['avg_abc'] = bin_score_df.groupby('bin_num')['abc_score'].transform('mean')
    bin_score_df_unique=bin_score_df.drop_duplicates(subset=('bin_num', 'avg_abc'))
    
    #initiate weight track of zeros (no enhancer, not included in score )
    background_weight = 0
    weight_track = np.array([background_weight]*len(disruption_track), dtype=float)
    
    #loop through bins, and abc score to weights track
    for row in bin_score_df_unique.itertuples():
        weight_track[row.bin_num] = row.avg_abc
    #print('weight track', weight_track[25:30])
        
    # Add nans to weight track
    weight_track[np.isnan(disruption_track)] = 0

    if len(roi_in_map) == 0:
        roi_id_values = ['']
    else:
        roi_id_values = roi_in_map.roi_id.values

    # Mask nan values in disruption track
    disruption_track_masked = np.ma.MaskedArray(disruption_track, mask=np.isnan(disruption_track))

    return np.ma.average(disruption_track_masked, weights=weight_track), weight_track, roi_id_values



def get_genes_in_map(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN, genome = 'hg38'):
    
    
    '''
    Get protein coding genes in map region.
    
    '''

    # For duplications, we want to plot genes that correspond to the masked REF matrix
    if SVTYPE == 'DUP':
        map_start_coord = map_start_coord[1]
        rel_pos_map = rel_pos_map[0]
    else:
        map_start_coord = map_start_coord[0]
        rel_pos_map = rel_pos_map[1]
        
    map_region_BED = BedTool.from_dataframe(pd.DataFrame({'CHR' : [CHR],
                                                          'start' : [map_start_coord],
                                                          'end' : [map_start_coord + map_length]}))

    if genome == 'hg38':
        gene_annot_BED = gene_annot_hg38_BED
    elif genome == 'hg19':
        gene_annot_BED = gene_annot_hg19_BED
        
    
    genes_in_map_BED = map_region_BED.intersect(gene_annot_BED, wa = True, wb = True)

    if genes_in_map_BED == '':
        genes_in_map = pd.DataFrame(columns=['Start', 'End', 'Gene', 'width'])

    else:
        genes_in_map = (genes_in_map_BED
                        .to_dataframe()
                        .rename(columns = {'score': 'Start', 'strand':'End', 'thickStart':'Gene'})
                        [['Start', 'End', 'Gene']])
        
        genes_in_map.Start = [math.ceil(x/bin_size) for x in (genes_in_map.Start - map_start_coord)]
        genes_in_map.loc[genes_in_map.Start < 0,'Start'] = 0
        genes_in_map = genes_in_map[genes_in_map.Start <= 447]
        genes_in_map.End = [math.ceil(x/bin_size) for x in (genes_in_map.End - map_start_coord)]
        genes_in_map.loc[genes_in_map.End > target_length_cropped-1,'End'] = target_length_cropped-1
        
        if SVTYPE == 'DUP' and abs(int(SVLEN)) > bin_size/2:
            
            genes_in_map.loc[(genes_in_map.Start >= rel_pos_map) & 
                              (genes_in_map.End > rel_pos_map),
                              ['Start', 'End']] = genes_in_map.loc[(genes_in_map.Start > rel_pos_map) & 
                                                                  (genes_in_map.End > rel_pos_map),
                                                                  ['Start', 'End']] + math.ceil(SVLEN/bin_size)


            middle_genes = genes_in_map.loc[(genes_in_map.Start < rel_pos_map) & 
                              (genes_in_map.End > rel_pos_map)].index


            for middle_gene in middle_genes:

                left = genes_in_map.loc[[middle_gene]]
                left.End = rel_pos_map
                right = genes_in_map.loc[[middle_gene]]
                right.Start = rel_pos_map + math.ceil(SVLEN/bin_size)
                right.End = right.End + math.ceil(SVLEN/bin_size)

                genes_in_map = pd.concat([genes_in_map.drop(middle_gene), left, right], axis = 0).reset_index(drop=True)
            
            # Remove genes that left prediction window
            genes_in_map = genes_in_map[~((genes_in_map.End > target_length_cropped) &
                                        (genes_in_map.Start > target_length_cropped))]
            
            # Crop end of genes that left prediction window
            genes_in_map.loc[(genes_in_map.End > target_length_cropped),'End'] = target_length_cropped
            
            genes_in_map.reset_index(drop=True, inplace = True)

    
        genes_in_map['width'] = [x-y+1 for x,y in zip(genes_in_map.End, genes_in_map.Start)]
            
            
    return genes_in_map



def plot_maps_abc(maps, genes_in_map, lines, roi_in_map, disruption_track, weight_track, title, scale = 1, gene_rows = 3):
    
    '''
    Plot the reference and alternate predicted contact frequency maps with lines at the beginning and end of the variant. 
    For duplications, there will be 3 lines marking the two regions that are duplicates.    
    Plot genes that match the regions in the reference map.
    Plot disruption score tracks.
    
    '''

    if len(genes_in_map)>0:
        gene_track = list(genes_in_map[['Start', 'width']].to_records(index = False))
    roi_track = list(roi_in_map[['Start', 'width']].to_records(index = False))
    
    n = maps[0].shape[0]

    
    linestyle = 'dashed'
    linewidth = 1.5*scale
    plot_width = 2*scale
    plot_width1D = 1*scale
    #plot_width1D determines height of the disruption, activity, etc tracks 
    fig, gs = gridspec_inches([plot_width*2], [plot_width, .15, plot_width, .15, plot_width, .15, 
                                               plot_width1D, .15, plot_width1D, .15, plot_width1D, .15, plot_width1D/2])


    for i in range(len(maps)):
        
        plt.subplot(gs[i*2,0])
        #map_plot=mat_from_vector(maps[i])
        #plot_upright(plt, map_plot, lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)
        plotting_utils.pcolormesh_45deg(plt, maps[i], lines, linewidth, cmap= 'RdBu_r', vmax=2, vmin=-2)

     
    #gets difference map 
    #print('plotting difference ')
    plt.subplot(gs[4,0])
    plotting_utils.pcolormesh_45deg(plt, maps[0] - maps[1], lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)
    #map_diff_plot=mat_from_vector(maps[0] - maps[1])

    #plot_upright(plt, map_diff_plot, lines, linewidth, cmap= 'PRGn', vmax=1, vmin=-1)



    # Disruption track
   # print ('plotting disruption')
    plt.subplot(gs[6,0])
    plt.plot(list(range(target_length_cropped)), disruption_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.ylabel('Disruption') 
    #plt.yticks([])


#     # ROI regions track
#     ax1 = fig.add_subplot(gs[8,0])
            
#     ax1.broken_barh(roi_track, (0, 1*scale), facecolors='tab:gray')
       
#     if lines is not None:
#         for line in lines:
#             plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
#     plt.ylim([0,1])
#     plt.xlim([0,448])
#     plt.xticks([])
#     plt.yticks([])


    # Weight track
    #print('plotting weights')
    plt.subplot(gs[8,0])
    plt.plot(list(range(target_length_cropped)), weight_track, color = 'black', linewidth = linewidth)
    plt.xlim([0,target_length_cropped])
    plt.gca().yaxis.tick_right()
    if lines is not None:
        for line in lines:
            plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
    plt.xticks([])
    plt.ylabel('ABC Activity') 

    #plt.yticks([])
    

    # Gene track
    if len(genes_in_map)>0:
        ax1 = fig.add_subplot(gs[10,0])


        gene_bar_height = 15*gene_rows/3
        bar_locations = []

        for bar_level in range(gene_rows):

            bar_location = (bar_level+1)*(gene_bar_height/gene_rows) - 0.5
            ax1.broken_barh(gene_track[bar_level:][::gene_rows], 
                            (bar_location, 0.5*scale), 
                            facecolors='tab:blue')
            bar_locations.append(bar_location)


        bar_locations = (bar_locations*math.ceil(len(genes_in_map)/gene_rows))[:len(genes_in_map)]

        for i in range(len(genes_in_map)):

            bar_location = bar_locations[i]

            gene = genes_in_map.loc[i,'Gene']
            location = genes_in_map.loc[i,'Start']

            ax1.annotate(gene, (location,bar_location), # annotate gene at the start
                         rotation = 45, 
                         ha = "right", va = "top", # set horizontal and vertical alignment
                         annotation_clip = False, # keep annotation if outside of window
                         fontsize = 8*scale) 
            
        if lines is not None:
            for line in lines:
                plt.axvline(x=line, color='gray', linestyle=linestyle, linewidth = linewidth)
        plt.ylabel('Genes',rotation=90)
        plt.ylim([0,gene_bar_height])
        plt.xlim([0,448])
        plt.axis('off')
        
    
    fig.suptitle(title, fontsize=12, y=1.12)


    #plt.show()
    return(fig)


#######################################################
#RUN ABC WEIGHTED SCORING FOR EACH SET OF SAMPLES 



def run_abc_scores(scores, full_scores, maps, tracks, abc_BED):
    # Plot maps for the chosen variant
    figures=[]

    for var_index in range(0, len(scores)):
    #for var_index in range(0,3):

        try:

            print('var index', var_index)
            variant_id = f'{var_index}_maps_HFF_0'


            REF_pred, ALT_pred = maps[variant_id][:2]
            REF_pred, ALT_pred = [get_Akita_scores_utils.map_from_vector(x) for x in [REF_pred, ALT_pred]]

            #rel_pos_map, map_start_coord = maps[variant_id][2:]
            #print(maps[variant_id])
            rel_pos_map = maps[variant_id][2]
            map_start_coord = maps[variant_id][3]

            print('rel pos map', rel_pos_map)
            #print('map start coord', map_start_coord)
            disruption_mse_track=tracks[f'{var_index}_mse_track_HFF_0']
            #print(rel_pos_map)


            CHR = scores.iloc[var_index].CHROM
            start = scores.iloc[var_index].POS
            end = scores.iloc[var_index].END
            SVTYPE = scores.iloc[var_index].SVTYPE



            if 'SVLEN' in scores.columns:
                SVLEN = scores.iloc[var_index].SVLEN
            else:
                SVLEN = '-'


            if SVTYPE=='BND':
                #print('running BND')

                REF=scores.iloc[var_index].REF
                ALT=scores.iloc[var_index].ALT
                #print('coords', (CHR, start, end))
                #print('ALT is ', ALT)


                #EDIT -- CALL BND COORDINATES WITHIN EACH FUNCTION 
                #KEEP GET ROI IN MAP CLEANER        

                #this gets the lines for where the SV is 
                lines = plotting_utils.get_var_bins(rel_pos_map, SVTYPE, SVLEN)

                #gets which enhancer regions are in the map 
                roi_in_map= get_roi_in_map_BND((CHR, start, REF, ALT, 0, False), rel_pos_map[0], abc_BED)
                #print('num enahncer in map', len(roi_in_map))

                #gets weighted ABC score, using A from ABC model, C from Akita prediction 
                weighted_scores=get_weighted_ABC_score_BND(disruption_mse_track, roi_in_map)
                weighted_abc_score=weighted_scores[0]
                #print('weighted_scores', weighted_abc_score)


                #add weighted scores, and number of enhancer regions in map to dataframe 
                orig_index=scores.index[var_index]
                full_scores.loc[orig_index, 'ABC_score'] = weighted_abc_score
                mse_score=full_scores.loc[orig_index].mse_mean
                #print(mse_score)

                #get the genes in the map region
                #genes_in_map=plotting_utils.get_genes_in_map_BND(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN)
                genes_in_map = get_genes_in_map_BND((CHR, start, REF, ALT, 0, False), rel_pos_map[0])


                BND_left_coord, BND_right_coord, BND_rel_pos, antisense = get_coords_BND(CHR, start, REF, ALT, 0, False)

                #title for figure
                title=f"""
                Variant {orig_index}, {SVTYPE} 
                BND Left: {BND_left_coord[0]}:{(BND_left_coord[1]):,}-{(BND_left_coord[2]):,}
                BND Right: {BND_right_coord[0]}:{(BND_right_coord[1]):,}-{(BND_right_coord[2]):,}
                MSE Disruption: {mse_score:.5f} 
                Weighted ABC: {weighted_abc_score:.5f} 
                # Enhancers in Map: {len(roi_in_map)}"""

                #plot
                fig_var=plot_maps_abc_BND([REF_pred, ALT_pred], genes_in_map, lines, roi_in_map, 
                                      disruption_mse_track, weighted_scores[1], title=title, scale = 1, gene_rows = 4)
                figures.append(fig_var)

            else: #SVTYPE is DEL, DUP, INV   
                #this gets the lines for where the SV is 
                lines = plotting_utils.get_var_bins(rel_pos_map, SVTYPE, SVLEN)

                #gets which enhancer regions are in the map 
                roi_in_map=get_roi_in_map(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN, abc_BED)
                #print('num enahncer in map', len(roi_in_map))

                #gets weighted ABC score, using A from ABC model, C from Akita prediction 
                weighted_scores=get_weighted_ABC_score(disruption_mse_track, roi_in_map)
                weighted_abc_score=weighted_scores[0]
                #print('weighted_scores', weighted_abc_score)


                #add weighted scores, and number of enhancer regions in map to dataframe 
                orig_index=scores.index[var_index]
                full_scores.loc[orig_index, 'ABC_score'] = weighted_abc_score
                mse_score=full_scores.loc[orig_index].mse_mean
                #print(mse_score)

                #get the genes in the map region
                genes_in_map=get_genes_in_map(CHR, map_start_coord, rel_pos_map, SVTYPE, SVLEN)


                #title for figure
                title=f'Variant {orig_index}\n{SVTYPE} at {CHR}:{(start):,}-{(end):,}\n MSE Disruption: {mse_score:.5f} \n Weighted ABC: {weighted_abc_score:.5f} \n # Enhancers in Map: {len(roi_in_map)}'

                #plot
                fig_var=plot_maps_abc([REF_pred, ALT_pred], genes_in_map, lines, roi_in_map, 
                                      disruption_mse_track, weighted_scores[1], title=title, scale = 1, gene_rows = 4)
                figures.append(fig_var)

        except:
            print('error in var index', var_index)
            
    return(figures, full_scores)

        

#######################################################
#READ IN ABC AND SCORES DATA FOR EACH TUMOR TYPE  


out_dir='/pollard/data/projects/shzhang/cbtn/results/svs_tumortype_abc_FINAL/'
abc_res_dir='/pollard/home/shzhang/ABC-Enhancer-Gene-Prediction/results/'

tt_list=['atrt', 'dipg', 'medullo', 'phgg']
tt_abc_list=['BT16', 'DIPG007_2', 'D283', 'KNS42']


for a in range(0,4):
    tt=tt_list[a]
    tt_abc=tt_abc_list[a]

    #read in abc scores
    abc_scores_tt=pd.read_csv(f'{abc_res_dir}{tt_abc}/Predictions/EnhancerPredictionsFull_threshold0.025_self_promoter.tsv', sep='\t')


    #make an abc roi bedtool file
    roi_file=abc_scores_tt
    roi_file['roi_id'] = ['roi_' + str(x) for x in roi_file.index]
    single_cols = ['chr', 'start', 'end', 'roi_id', 'activity_base']
    roi_coords = roi_file[single_cols]
    abc_tt_BED=BedTool.from_dataframe(roi_coords.astype({'start': int, 'end': int}))



    tt_file=f"{tt}_top10mse"

    maps_tt = np.load(f'{out_dir}{tt}/{tt_file}_maps.npy', allow_pickle="TRUE").item()
    tracks_tt = np.load(f'{out_dir}{tt}/{tt_file}_tracks.npy', allow_pickle="TRUE").item()

    scores_tt=pd.read_csv(f'{out_dir}{tt}/{tt_file}.txt', sep='\t', index_col=0)
    scores_tt['var_index']=range(len(scores_tt))

    full_scores_tt=CBTN_SV_scores.loc[scores_tt.index]
    print(f'number of {tt} SVs: {len(full_scores_tt)}')
    
    #actually run the abc scores 
    figures,full_scores= run_abc_scores(scores_tt, full_scores_tt, maps_tt, tracks_tt, abc_tt_BED)

    
    #calculate number of files needed to save figures in increments of 500
    n_fig_files=math.ceil(len(full_scores_tt)/500)
    print('n fig filer', n_fig_files)
    
    #save figures in increments 1000, otherwise big 
    for index in range(0, n_fig_files):
        i_start=index*500
        figures_subset=figures[i_start:(i_start+500)]
        save_figures_to_pdf(f'{out_dir}{tt}/abc/{tt_file}_abc_plots_{i_start}.pdf', figures_subset)

    
    #save scores df with new abc scores 
    full_scores.to_csv(f'{out_dir}{tt}/abc/{tt_file}_abc_scores.csv')




