import os
import pandas as pd
from operator import xor

# unused
# from numpy import mean
# mean([s1[i]==s2[i] for i in xrange(len(s1))])
seq1 = [char for char in "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG-------------------"]
seq2 = [char for char in "------CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG-------------"]
seq3 = [char for char in "-------------------ACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGCCATG"]

MSA = {"sequence_1":seq1, "sequence_2":seq2, "sequence_3": seq3}
MSA_df = pd.DataFrame()

for key in MSA:
    if key != 'alignment':
        key_index = []
        count = 1
        for letter in MSA[key]:
            if letter != "-":
                key_index.append(count)
                count += 1
            else:
                key_index.append(' ')
        MSA_df[key+'_index'] = key_index
    MSA_df[key] = MSA[key]
print(MSA_df)
print(MSA_df.columns)

print(MSA_df['sequence_2_index'])
os.chdir('c:/Users/fejsa/OneDrive/Desktop/Graduation project/COMSAA/input')
CDS = pd.read_csv('test_coding.csv')
CDS_start = CDS['Start']
CDS_stop = CDS['Stop']
CDS_strand = CDS['Strand']

# Separate the sequence into coding and non-coding, so the difference between them could be checked
CDS_combined = set()

for i in CDS_start.index:
    CDS_combined = CDS_combined.union(
            set(range(CDS_start[i], CDS_stop[i]+1)))
print(CDS_combined)
CDS_combined_df = MSA_df[MSA_df.iloc[:, 0].isin(CDS_combined)]
nonCDS_combined_df = MSA_df[MSA_df.iloc[:, 0].isin(
        CDS_combined) == False]
CDS_final = CDS_combined_df.iloc[:, range(1, CDS_combined_df.shape[1],2)]
nonCDS_final = nonCDS_combined_df.iloc[:, range(1, nonCDS_combined_df.shape[1],2)]
names = ['sequence_1', 'sequence_2', 'sequence_3']
'''
cds combined df - index 456789 10 11 12 16 17 18 19 20 21 22 23 .. 44? 
treba do 57?
noncds combined df
ostalo

'''
# make new dataframes which will contain different calculations types
# percent identity matrix
nonCDS_similarity = pd.DataFrame(columns=names, index=names)
# mutation number matrix
nonCDS_mutation_num = pd.DataFrame(columns=names, index=names)
# transitions number matrix
nonCDS_transition_num = pd.DataFrame(columns=names, index=names)
# transversions number matrix
nonCDS_transversion_num = pd.DataFrame(columns=names, index=names)
# transition/transvertion ratio matrix
nonCDS_tt_ratio = pd.DataFrame(columns=names, index=names)
# number of gaps matrix
nonCDS_gaps = pd.DataFrame(columns=names, index=names)
# number of indels matrices
nonCDS_insertions = pd.DataFrame(columns=names, index=names)
nonCDS_deletions = pd.DataFrame(columns=names, index=names)

for name1 in names:
    seq1 = nonCDS_final[name1].values
    for name2 in names:
        seq2 = nonCDS_final[name2].values
        length = len(seq1)
        # initialize values to be added
        mutations = 0
        transitions = 0
        transversions = 0
        gaps = 0
        insertions = 0
        deletions = 0
        Ns = 0
        # track every nucleotide from both sequences
        for i in range(length):
            nuc1 = seq1[i]
            nuc2 = seq2[i]
            # only take into account the differences
            if nuc1 == 'N' or nuc2 == 'N':
                Ns += 1
            elif nuc1 != nuc2:
                mutations += 1
                if nuc1 == '-':
                    gaps += 1
                    insertions += 1
                elif nuc2 == '-':
                    gaps += 1
                    deletions += 1
                elif (nuc1 == 'A' and nuc2 == 'G') or (nuc1 == 'G' and nuc2 == 'A'):
                    transitions += 1
                elif (nuc1 == 'T' and nuc2 == 'C') or (nuc1 == 'C' and nuc2 == 'T'):
                    transitions += 1
                else:
                    transversions += 1

        # populate matrices
        nonCDS_similarity.at[name1, name2] = 1-mutations/length
        nonCDS_mutation_num.at[name1, name2] = mutations
        nonCDS_transition_num.at[name1, name2] = transitions
        nonCDS_transversion_num.at[name1, name2] = transversions
        try:
            nonCDS_tt_ratio.at[name1,
                               name2] = transitions/transversions
        except:
            nonCDS_tt_ratio.at[name1, name2] = 0
        nonCDS_gaps.at[name1, name2] = gaps
        nonCDS_insertions.at[name1, name2] = insertions
        nonCDS_deletions.at[name1, name2] = deletions


print(nonCDS_insertions)
print(nonCDS_deletions)
