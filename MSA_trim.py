#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd

#read in data
df = pd.read_csv('test.csv')


# In[3]:


ref_list = df.loc[0][1:df.shape[1]].tolist()


# In[4]:


# scan through each column, ignore the first column
snp_columns = list(df)
# let's get information about the reference sequence too

#i is an index 
#extract no. of samples
column_num = 0
sample_num = 0
ambig_num = 0
ref_num = 0
indel_num = 0
snp_num = 0
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
                
        if len(snp_prop) == 1:
            del df[i]
        
        elif len(snp_prop) == 2:
            if (ambig_num/sample_num > 0.5) or (indel_num/sample_num > 0.5) or (snp_num/sample_num > 0.9):
                del df[i]
            else:
                continue
                
        else:
            if ((ambig_num + indel_num)/sample_num > 0.5) or (indel_num/sample_num > 0.5) or (ambig_num/sample_num > 0.5):
                del df[i]
            else:
                continue
    ambig_num = 0
    ref_num = 0
    indel_num = 0
    snp_num = 0        
    snp_prop.clear()
    column_num += 1
        


# In[39]:


#df.drop(df.index[[len(df[i])-1,len(df[i])-2]])
for index,row in df.iterrows():
    if row["X.reference_pos"] == "map-quality" or row["X.reference_pos"] == "annotations":
        continue
    else:
        print(">{}".format(row["X.reference_pos"]))
        seq = list(row[1:len(row)].to_string(index=False))
        seq = [e for e in seq if e not in (' ', '\n')]
        str_seq = ''.join(map(str, seq))
        print(str_seq)

