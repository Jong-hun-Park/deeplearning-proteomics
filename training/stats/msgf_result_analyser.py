import os
import sys

import numpy as np
import scipy.stats as stats
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from astropy.units import percent

'''
Created on 2017. 1. 10.

@author: JONGHUN
'''

# USER parameter
RESULT_EXTENSION = '.tsv'  # result file name extension
DECOY_LABEL = 'XXX_'

target_filenames = []

def search_files(dirname):
    filenames = os.listdir(dirname)
       
    for filename in filenames:
        full_filename = os.path.join(dirname, filename)
        
        if os.path.isdir(full_filename):
            search_files(full_filename)
             
        ext = os.path.splitext(full_filename)[-1]
        
        # target files only having specific extension
        if ext == RESULT_EXTENSION:
            target_filenames.append(full_filename)
            
    return target_filenames
            
def read_msgf_result(input_file):
    result = Msgf_result()
    
    # For header
    line = input_file.readline()
    
    while (True):
        line = input_file.readline()
        if not line: break
        
        psm = Msgf_psm(line)
        result.add_psm(psm)
    
    return result
class Msgf_psm():
    
    # TODO: the index can be varied by the search option. so, should be parameterized 
    def __init__(self, line):
        splited = line.split("\t")
        
        if (len(splited) > 14):
            self.specFile = splited[0]
            self.specID = splited[1]
            self.scanNum = splited[2]
            self.title = splited[3]
            self.fragMethod = splited[4]
            self.precursor = splited[5]
            self.isotopeError = splited[6]
            self.precursorError = splited[7]
            self.charge = splited[8]
            self.peptide = splited[9]
            self.protein = splited[10]
            self.denovoScore = splited[11]
            self.msgfScore = splited[12]
            self.specEvalue = splited[13]
            self.EValue = splited[14]
            self.QValue = splited[15]
            self.pepQValue = splited[16]
            
            self.is_unmodi = Msgf_psm.isunmodi(splited[9])
            self.strip_peptide = Msgf_psm.get_strip_peptide(splited[9])
        
        else:
            print "this result has not applied FDR (no Q-value column)"
            SystemExit
        
        def isunmodi(peptide):
            return not any(char.isdigit() for char in peptide)
        
        def get_strip_peptide(peptide):
            strip_peptide = ""
            for aa in peptide:
                if aa.isalpha():
                    strip_peptide += aa 
            
            return strip_peptide
            
class Msgf_result():
# A class for each result file.
    PSM_count = 0
    
    target_count = 0
    decoy_count = 0
    
    charge_2_count = 0
    charge_3_count = 0
    charge_more_than_4_count = 0
    
    pep_len_6_count = 0
    pep_len_7_count = 0
    pep_len_8_count = 0
    pep_len_9_count = 0
    pep_len_10_count = 0
    pep_len_11_count = 0
    pep_len_12_count = 0
    pep_len_13_count = 0
    pep_len_14_count = 0
    pep_len_15_count = 0
    pep_len_more_than_16_count = 0
    
    
    def __init__(self, qvalue_cutoff):
        self.PSM_count = 0
        self.unmodi_count = 0
        
        self.target_count = 0
        self.decoy_count = 0
        
        # Charge state
        self.charge_2_count = 0
        self.charge_3_count = 0
        self.charge_more_than_4_count = 0
        
        self.pep_len_6_count = 0
        self.pep_len_7_count = 0
        self.pep_len_8_count = 0
        self.pep_len_9_count = 0
        self.pep_len_10_count = 0
        self.pep_len_11_count = 0
        self.pep_len_12_count = 0
        self.pep_len_13_count = 0
        self.pep_len_14_count = 0
        self.pep_len_15_count = 0
        self.pep_len_more_than_16_count = 0
        
        self.qvalue_cutoff = qvalue_cutoff
        
        self.psmlist = []
        self.scorelist = []
                
# TODO: plotting function
# TODO: Q-value cutoff
    def get_stats(self, td, cs, length, score, qvalue_cutoff):
        
        for psm in self.psmlist:
            if (psm.QValue <= self.qvalue_cutoff):
                if td:
                    if (psm.protein.startswith(DECOY_LABEL)):
                        self.decoy_count += 1
                        Msgf_result.decoy_count += 1
                    else:
                        self.target_count += 1
                        Msgf_result.target_count += 1
                
                    # Error check
                    td_sum = self.target_count + self.decoy_count
                    if (td_sum != self.PSM_count):
                        print ("TD error")
        
                if cs:
                    # draw charge distribution graph
                    if (psm.charge == '2'):
                        self.charge_2_count += 1
                        Msgf_result.charge_2_count += 1
                    elif (psm.charge == '3'):
                        self.charge_3_count += 1
                        Msgf_result.charge_3_count += 1
                    else:
                        self.charge_more_than_4_count += 1
                        Msgf_result.charge_more_than_4_count += 1
            
                if length:
                    # draw length distribution graph
                    if (len(psm.strip_peptide) == 6):
                        self.pep_len_6_count += 1
                        Msgf_result.pep_len_6_count += 1
                    elif (len(psm.strip_peptide) == 7):
                        self.pep_len_7_count += 1
                        Msgf_result.pep_len_7_count += 1
                    elif (len(psm.strip_peptide) == 8):
                        self.pep_len_8_count += 1
                        Msgf_result.pep_len_8_count += 1
                    elif (len(psm.strip_peptide) == 9):
                        self.pep_len_9_count += 1
                        Msgf_result.pep_len_9_count += 1
                    elif (len(psm.strip_peptide) == 10):
                        self.pep_len_10_count += 1
                        Msgf_result.pep_len_10_count += 1 
                    elif (len(psm.strip_peptide) == 11):
                        self.pep_len_11_count += 1
                        Msgf_result.pep_len_11_count += 1
                    elif (len(psm.strip_peptide) == 12):
                        self.pep_len_12_count += 1
                        Msgf_result.pep_len_12_count += 1
                    elif (len(psm.strip_peptide) == 13):
                        self.pep_len_13_count += 1
                        Msgf_result.pep_len_13_count += 1
                    elif (len(psm.strip_peptide) == 14):
                        self.pep_len_14_count += 1
                        Msgf_result.pep_len_14_count += 1
                    elif (len(psm.strip_peptide) == 15):
                        self.pep_len_15_count += 1
                        Msgf_result.pep_len_15_count += 1
                    else:
                        self.pep_len_more_than_16_count += 1
                        Msgf_result.pep_len_more_than_16_count += 1
                
                if score:
        #             # draw a score distribution graph
        #             self.scorelist.append(float(psm.EValue))
        #             
        #             h = sorted(self.scorelist)
        #             fit = stats.norm.pdf(h, np.mean(h), np.std(h))  #this is a fitting indeed
        # 
        #             plt.plot(h,fit,'-o')
        #             
        #             plt.hist(h, bins=10, normed=True)      #use this to draw histogram of your data
        #             
        #             plt.show()  
                    pass
            
    def print_stats(self, td, cs, length, score):
        print "-----PSM count-----"
        print self.PSM_count
        print "unmodi count:\t" + str(self.unmodi_count)
        
        if td:
            print "-----Target and decoy count-----"
            print "target count:\t"  + str(self.target_count)
            print "decoy count:\t"  + str(self.decoy_count)
        
        if cs:
            print "-----Charge state-----"
            print "charge 2+:\t" + str(self.charge_2_count)
            print "charge 3+:\t" + str(self.charge_3_count)
            print "charge >= 4+:\t" + str(self.charge_more_than_4_count)
        
        if length:
            print "-----Peptide length-----"
            print "length 6:\t" + str(self.pep_len_6_count)
            print "length 7:\t" + str(self.pep_len_7_count)
            print "length 8:\t" + str(self.pep_len_8_count)
            print "length 9:\t" + str(self.pep_len_9_count)
            print "length 10:\t" + str(self.pep_len_10_count)
            print "length 11:\t" + str(self.pep_len_11_count)
            print "length 12:\t" + str(self.pep_len_12_count)
            print "length 13:\t" + str(self.pep_len_13_count)
            print "length 14:\t" + str(self.pep_len_14_count)
            print "length 15:\t" + str(self.pep_len_15_count)
            print "length >= 16:\t" + str(self.pep_len_more_than_16_count)
            
        if score:
            pass
        
        print "\n"
    
    def write_stats(self, output_file, td, cs, length, score):
        output_file.write("-----PSM count-----" + '\n')
        output_file.write(str(self.PSM_count) + '\n')
        output_file.write("unmodi count:\t" + str(self.unmodi_count) + '\n')
        
        if td:
            output_file.write("-----Target and decoy count-----" + '\n')
            output_file.write("target count:\t"  + str(self.target_count) + '\n')
            output_file.write("decoy count:\t"  + str(self.decoy_count) + '\n')
        
        if cs:
            output_file.write("-----Charge state-----" + '\n')
            output_file.write("charge 2+:\t" + str(self.charge_2_count) + '\n')
            output_file.write("charge 3+:\t" + str(self.charge_3_count) + '\n')
            output_file.write("charge >= 4+:\t" + str(self.charge_more_than_4_count) + '\n')
        
        if length:
            output_file.write("-----Peptide length-----" + '\n')
            output_file.write("length 6:\t" + str(self.pep_len_6_count) + '\n')
            output_file.write("length 7:\t" + str(self.pep_len_7_count) + '\n')
            output_file.write("length 8:\t" + str(self.pep_len_8_count) + '\n')
            output_file.write("length 9:\t" + str(self.pep_len_9_count) + '\n')
            output_file.write("length 10:\t" + str(self.pep_len_10_count) + '\n')
            output_file.write("length 11:\t" + str(self.pep_len_11_count) + '\n')
            output_file.write("length 12:\t" + str(self.pep_len_12_count) + '\n')
            output_file.write("length 13:\t" + str(self.pep_len_13_count) + '\n')
            output_file.write("length 14:\t" + str(self.pep_len_14_count) + '\n')
            output_file.write("length 15:\t" + str(self.pep_len_15_count) + '\n')
            output_file.write("length >= 16:\t" + str(self.pep_len_more_than_16_count) + '\n')
            
        if score:
            pass
    
    def print_total_stats(self, td, cs, length, score):
        print target_filenames
        print "-----PSM count-----"
        print Msgf_result.PSM_count
        print Msgf_result.unmodi_count
        
        if td:
            print "-----Target and decoy count-----"
            print "target count: "  + str(Msgf_result.target_count)
            print "decoy count: "  + str(Msgf_result.decoy_count)
        
        if cs:
            print "-----Charge state-----"
            print "charge 2+:    " + str(Msgf_result.charge_2_count)
            print "charge 3+:    " + str(Msgf_result.charge_3_count)
            print "charge >= 4+: " + str(Msgf_result.charge_more_than_4_count)
        
        if length:
            print "-----Peptide length-----"
            print "length 6:     " + str(Msgf_result.pep_len_6_count)
            print "length 7:     " + str(Msgf_result.pep_len_7_count)
            print "length 8:     " + str(Msgf_result.pep_len_8_count)
            print "length 9:     " + str(Msgf_result.pep_len_9_count)
            print "length 10:    " + str(Msgf_result.pep_len_10_count)
            print "length 11:    " + str(Msgf_result.pep_len_11_count)
            print "length 12:    " + str(Msgf_result.pep_len_12_count)
            print "length 13:    " + str(Msgf_result.pep_len_13_count)
            print "length 14:    " + str(Msgf_result.pep_len_14_count)
            print "length 15:    " + str(Msgf_result.pep_len_15_count)
            print "length >= 16: " + str(Msgf_result.pep_len_more_than_16_count)
            
        if score:
            pass
        
        print "\n"
        
    def write_total_stats(self, output_file, td, cs, length, score):
            
        output_file.write(target_filenames)
        output_file.write('\n')
        
        output_file.write("-----PSM count-----" + '\n')
        output_file.write(str(Msgf_result.PSM_count) + '\n')
        output_file.write("unmodi count:\t" + str(Msgf_result.unmodi_count) + '\n')
        
        if td:
            output_file.write("-----Target and decoy count-----" + '\n')
            output_file.write("target count:\t"  + str(Msgf_result.target_count) + '\n')
            output_file.write("decoy count:\t"  + str(Msgf_result.decoy_count) + '\n')
        
        if cs:
            output_file.write("-----Charge state-----" + '\n')
            output_file.write("charge 2+:\t" + str(Msgf_result.charge_2_count) + '\n')
            output_file.write("charge 3+:\t" + str(Msgf_result.charge_3_count) + '\n')
            output_file.write("charge >= 4+:\t" + str(Msgf_result.charge_more_than_4_count) + '\n')
        
        if length:
            output_file.write("-----Peptide length-----" + '\n')
            output_file.write("length 6:\t" + str(Msgf_result.pep_len_6_count) + '\n')
            output_file.write("length 7:\t" + str(Msgf_result.pep_len_7_count) + '\n')
            output_file.write("length 8:\t" + str(Msgf_result.pep_len_8_count) + '\n')
            output_file.write("length 9:\t" + str(Msgf_result.pep_len_9_count) + '\n')
            output_file.write("length 10:\t" + str(Msgf_result.pep_len_10_count) + '\n')
            output_file.write("length 11:\t" + str(Msgf_result.pep_len_11_count) + '\n')
            output_file.write("length 12:\t" + str(Msgf_result.pep_len_12_count) + '\n')
            output_file.write("length 13:\t" + str(Msgf_result.pep_len_13_count) + '\n')
            output_file.write("length 14:\t" + str(Msgf_result.pep_len_14_count) + '\n')
            output_file.write("length 15:\t" + str(Msgf_result.pep_len_15_count) + '\n')
            output_file.write("length >= 16:\t" + str(Msgf_result.pep_len_more_than_16_count) + '\n')
            
        if score:
            pass
            
    def add_psm(self, psm):
        self.psmlist.append(psm)
        
        self.PSM_count += 1
        Msgf_result.PSM_count += 1

        if (psm.is_unmodi):
            self.unmodi_count += 1

if __name__ == '__main__':
    qvalue_cutoff = 1
    
    selected_filenames = search_files("I:\CPTAC\mgf")
    
    if (len(selected_filenames) == 0):
        print "please check the result file extension"
        sys.exit(-1)
        
    for file_name in selected_filenames:
        print "File name : " + str(file_name)
        
        input_file = open(file_name, 'r')
        output_file = open(file_name + ".stats.out", 'w')
        
        msgf_result = read_msgf_result(input_file)
        msgf_result.get_stats(td=True, cs=True, length=True, score=True)
        msgf_result.print_stats(td=True, cs=True, length=True, score=False)
        msgf_result.write_stats(output_file, td=True, cs=True, length=True, score=False)
        
        input_file.close()
        output_file.close()
    
    total_output = open("stat_total.out", "w")
    msgf_result.print_total_stats(td=True, cs=True, length=True, score=False);
    msgf_result.write_total_stats(total_output, td=True, cs=True, length=True, score=False)
    print "Finished"
