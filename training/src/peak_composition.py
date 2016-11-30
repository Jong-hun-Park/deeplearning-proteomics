'''
Created on 2016. 11. 30.

@author: JONGHUN
'''

import aa_mass_table

aa_count = {
'A' : 0,
'C' : 0,
'D' : 0,
'E' : 0,
'F' : 0,
'G' : 0,
'H' : 0,
'I' : 0,
'K' : 0,
'L' : 0,
'M' : 0,
'N' : 0,
'P' : 0,
'Q' : 0,
'R' : 0,
'S' : 0,
'T' : 0,
'V' : 0,
'W' : 0,
'Y' : 0
}

'''
b ion - sum of AA from 1 to i-3 
y ion - sum of AA from i+3 to len(peptide)
'''
def get_nterm_feature(peptide, ion_type, ion_index):
    b_possible = True
    y_possible = True
    
    if ion_index < 3:
        b_possible = False
    if len(peptide) - ion_index < 3:
        y_possible = False
    
    if b_possible and ion_type == 'b':
        for i in range(ion_index - 3):
            aa_count[peptide[i]] += 1
            
    elif y_possible and ion_type == 'y':
        for i in range(ion_index + 3, len(peptide)):
            aa_count[peptide[i]] += 1
    
    return aa_count.values()
        
def get_cterm_feature(peptide, ion_type, ion_index):
    pass

def get_composition_feature(peptide, ion_type, ion_index):
    nterm_feature = get_nterm_feature(peptide, ion_type, ion_index)
    cterm_feature = get_cterm_feature(peptide, ion_type, ion_index)
    
    print nterm_feature



get_composition_feature("AAAAAAAAA", 'b', 5)