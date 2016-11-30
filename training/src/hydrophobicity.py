'''
Created on 2016. 11. 30.

@author: JONGHUN
'''
import hydrophobicity_table

def get_sum_hydf(peptide, ion_type, ion_index):
    sum_hydf = 0
    
    if ion_type == 'b':
        for i in range(ion_index):
            sum_hydf += hydrophobicity_table.get_aa_hydph(peptide[i])
    else:
        for i in range(ion_index, len(peptide)):
            sum_hydf += hydrophobicity_table.get_aa_hydph(peptide[i])
    
    return sum_hydf

def get_hydrphobicity_features(peptide, ion_type, ion_index):
    sum_hydf = get_sum_hydf(peptide, ion_type, ion_index)
    
    


        
     

