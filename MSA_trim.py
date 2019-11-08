#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import sys

#arguments
#1. the csv
csv = sys.argv[1]
#2. indel percentage
indel = float(sys.argv[2])
#3. ambiguity percentage
ambig = float(sys.argv[3])

#read in data
df = pd.read_csv(csv)
ref_list = df.loc[0][1:df.shape[1]].tolist()

# scan through each column, ignore the first column
snp_columns = list(df)
# let's get information about the reference sequence too

#i is an index 
#extract no. of samples
column_num = 0
sample_num = 0
ambig_num = 0
ambig_hist = []
ref_num = 0
ref_hist = []
indel_num = 0
indel_hist = []
snp_num = 0
snp_hist = []

snp_prop = {}
for i in snp_columns: 
    if "reference_pos" in i :
        continue
    else:  
        # able to get the column, now must change it to a list and 
        # cycle through by doing df[i][j].tolist() 
        snp_site = df[i][2:(len(df[i]))-2].tolist()
        # now we should create a python dictionary to associate SNPs and their proportions
        # amongst the samples.
        for snp in snp_site:
            if snp == "-":
                indel_num = indel_num +1
            elif snp != 'A' and snp != 'T' and snp != 'C' and snp != 'G':
                ambig_num = ambig_num + 1
            elif snp == ref_list[column_num]:
                ref_num = ref_num + 1
            else:
                snp_num = snp_num + 1
                
            if sample_num == 0:
                sample_num = len(snp_site)
                
            if snp_prop.get(snp) == None:
                snp_prop[snp] = 1
            else:
                snp_prop[snp] = snp_prop[snp]+1
        print("{},{},{},{}".format(indel_num,ambig_num,ref_num,snp_num))
        ambig_hist.append(ambig_num)
        ref_hist.append(ref_num)
        indel_hist.append(indel_num)
        snp_hist.append(snp_num) 
        
        if len(snp_prop) == 1:
            del df[i]
        
        elif len(snp_prop) == 2:
            if (ambig_num/sample_num > ambig) or (indel_num/sample_num > indel) or (snp_num/sample_num > 0.9):
                del df[i]
                
        else:
            if ((ambig_num + indel_num)/sample_num > 0.5) or (indel_num/sample_num > indel) or (ambig_num/sample_num > ambig):
                del df[i]
       
        ambig_num = 0
        ref_num = 0
        indel_num = 0
        snp_num = 0 
    
    snp_prop.clear()
    column_num += 1
        
