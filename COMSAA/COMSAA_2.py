# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 08:41:27 2021

@author: Muhamed
"""

import os
import pandas as pd
from operator import xor

# unused
# from numpy import mean
# mean([s1[i]==s2[i] for i in xrange(len(s1))])
os.chdir('c:/Users/fejsa/OneDrive/Desktop/Graduation project/COMSAA/')

for file in os.listdir('input'):
    print(file)
    if file.endswith('.clustal'):
        file_MSA = 'input/'+file
        file_CDS = file_MSA.replace('.clustal', '.csv')
        file_name = file.replace('.clustal', '')

        """
        Select input file
        """

        f = open(file_MSA, "r")

        # initiate a dictionary to be populated with sequences and alignment info from the file
        MSA = {}

        # skip the first three unnecessary lines of the input file
        for i in range(3):
            next(f)

        # determine at which position does the sequence start (as it can change from file to file, depending on the name length)
        first_line = f.readline().rstrip()
        index = first_line.index(first_line.split()[1])

        # return back to the start of the file and skip 3 unnecessary lines again
        f.seek(0)
        for i in range(3):
            next(f)

        # start reading every line
        for line in f:
            line = line.rstrip('\n')
            # skip empty lines
            if len(line) == 0:
                continue
            # copy alignment info over to the MSA dictionary (include lines with * and lines with only ' ')
            if '*' in line or len(line.strip()) == 0:
                if 'alignment' in MSA:
                    MSA['alignment'] += line[index:]
                else:
                    MSA['alignment'] = line[index:]
            # copy name and sequence info to the MSA dictionary
            else:
                name, sequence = line.split()[0], line.split()[1]
                if name in MSA:
                    MSA[name] += sequence
                else:
                    MSA[name] = sequence

        f.close()

        # change strings from the dictionary into lists of characters
        for key in MSA:
            MSA[key] = list(MSA[key])

        # make a dataframe which will contain all sequences together with their own respective indexes
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

        # make a dataframe which will only contain "mutated" segments
        MSA_mutations = MSA_df[MSA_df["alignment"]
                               == ' '].transpose().iloc[:-1]

        # generate the list of column/index names (except the 'alignment') for new data frames
        names = list(MSA.keys())
        names.remove('alignment')
        print('AAAA')
        print(MSA_df)
        # make new dataframes which will contain different calculations types
        # percent identity matrix
        MSA_similarity = pd.DataFrame(columns=names, index=names)
        # mutation number matrix
        MSA_mutation_num = pd.DataFrame(columns=names, index=names)
        # transitions number matrix
        MSA_transition_num = pd.DataFrame(columns=names, index=names)
        # transversions number matrix
        MSA_transversion_num = pd.DataFrame(columns=names, index=names)
        # transition/transvertion ratio matrix
        MSA_tt_ratio = pd.DataFrame(columns=names, index=names)
        # number of gaps matrix
        MSA_gaps = pd.DataFrame(columns=names, index=names)
        # number of indels matrices
        MSA_insertions = pd.DataFrame(columns=names, index=names)
        MSA_deletions = pd.DataFrame(columns=names, index=names)

        # TODO
        # only check mutated portions of the sequence to save on computational resources
        # but the length should stil refer to the original sequence

        # analyze all sequences with eachother
        for seq1 in MSA:
            # skip alignment sequences
            if seq1 == 'alignment':
                continue
            for seq2 in MSA:
                if seq2 == 'alignment':
                    continue
                length = len(MSA[seq1])
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
                    nuc1 = MSA[seq1][i]
                    nuc2 = MSA[seq2][i]
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
                    """ Uncomment 5 lines below if you want calculation to be done as Clustal Omega default """
                    # if nuc1 == '-' and nuc2 == '-':
                    #     length -= 1
                    # if xor(nuc1 == '-', nuc2 == '-'):
                    #     length -= 1
                    #     mutations -= 1
                    """ Uncomment 5 lines below if you want calculation to be done as Clustal Omega default """

                # populate matrices
                MSA_similarity.at[seq1, seq2] = 1-mutations/length
                MSA_mutation_num.at[seq1, seq2] = mutations
                MSA_transition_num.at[seq1, seq2] = transitions
                MSA_transversion_num.at[seq1, seq2] = transversions
                try:
                    MSA_tt_ratio.at[seq1, seq2] = transitions/transversions
                except:
                    MSA_tt_ratio.at[seq1, seq2] = 0
                MSA_gaps.at[seq1, seq2] = gaps
                MSA_insertions.at[seq1, seq2] = insertions
                MSA_deletions.at[seq1, seq2] = deletions

        print(MSA_deletions.to_string())

        # extract the info on the coding sequences from the second input file
        CDS = pd.read_csv(file_CDS)
        CDS_start = CDS['Start']
        CDS_stop = CDS['Stop']
        CDS_strand = CDS['Strand']

        # Separate the sequence into coding and non-coding, so the difference between them could be checked
        CDS_combined = set()

        for i in CDS_start.index:
            CDS_combined = CDS_combined.union(
                set(range(CDS_start[i], CDS_stop[i]+1)))
        print(len(CDS_combined))
        CDS_combined_df = MSA_df[MSA_df.iloc[:, 0].isin(CDS_combined)]
        nonCDS_combined_df = MSA_df[MSA_df.iloc[:, 0].isin(
            CDS_combined) == False]

        CDS_combined_df_mutations = CDS_combined_df[CDS_combined_df["alignment"]
                                                    == ' '].iloc[:, :-1]
        nonCDS_combined_df_mutations = nonCDS_combined_df[
            nonCDS_combined_df["alignment"] == ' '].iloc[:, :-1]

        CDS_final = CDS_combined_df.iloc[:, range(
            1, CDS_combined_df.shape[1], 2)]
        nonCDS_final = nonCDS_combined_df.iloc[:, range(
            1, nonCDS_combined_df.shape[1], 2)]

        # make new dataframes which will contain different calculations types
        # percent identity matrix
        CDS_similarity = pd.DataFrame(columns=names, index=names)
        # mutation number matrix
        CDS_mutation_num = pd.DataFrame(columns=names, index=names)
        # transitions number matrix
        CDS_transition_num = pd.DataFrame(columns=names, index=names)
        # transversions number matrix
        CDS_transversion_num = pd.DataFrame(columns=names, index=names)
        # transition/transvertion ratio matrix
        CDS_tt_ratio = pd.DataFrame(columns=names, index=names)
        # number of gaps matrix
        CDS_gaps = pd.DataFrame(columns=names, index=names)
        # number of indels matrices
        CDS_insertions = pd.DataFrame(columns=names, index=names)
        CDS_deletions = pd.DataFrame(columns=names, index=names)

        for name1 in names:
            seq1 = CDS_final[name1].values
            for name2 in names:
                seq2 = CDS_final[name2].values
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
                CDS_similarity.at[name1, name2] = 1-mutations/length
                CDS_mutation_num.at[name1, name2] = mutations
                CDS_transition_num.at[name1, name2] = transitions
                CDS_transversion_num.at[name1, name2] = transversions
                try:
                    CDS_tt_ratio.at[name1, name2] = transitions/transversions
                except:
                    CDS_tt_ratio.at[name1, name2] = 0
                CDS_gaps.at[name1, name2] = gaps
                CDS_insertions.at[name1, name2] = insertions
                CDS_deletions.at[name1, name2] = deletions

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

        print('NONCDS')

        print(type(nonCDS_final['sequence_1997'].values))
        for x in nonCDS_final['sequence_1997'].values:
            print(x)
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

        CDS_length = len(CDS_final)
        nonCDS_length = len(nonCDS_final)
        total_length = CDS_length + nonCDS_length
        CDS_nonCDS_ratio = CDS_length / nonCDS_length

        summary = pd.DataFrame(columns=['Similarity', 'Mutations', 'Mutation_Freq', 'Transitions', 'Transversions', 'TT_ratio', 'Gaps', 'Insertions', 'Deletions',
                                        'CDS_Similarity', 'CDS_Mutations', 'CDS_Mutation_Freq', 'CDS_Transitions', 'CDS_Transversions', 'CDS_TT_ratio', 'CDS_Gaps', 'CDS_Insertions', 'CDS_Deletions',
                                        'nonCDS_Similarity', 'nonCDS_Mutations', 'nonCDS_Mutation_Freq', 'nonCDS_Transitions', 'nonCDS_Transversions', 'nonCDS_TT_ratio', 'nonCDS_Gaps', 'nonCDS_Insertions', 'nonCDS_Deletions'])

        summary['Similarity'] = MSA_similarity.iloc[:, 0]
        summary['Mutations'] = MSA_mutation_num.iloc[:, 0]
        summary['Mutation_Freq'] = summary['Mutations']/total_length
        summary['Transitions'] = MSA_transition_num.iloc[:, 0]
        summary['Transversions'] = MSA_transversion_num.iloc[:, 0]
        summary['TT_ratio'] = MSA_tt_ratio.iloc[:, 0]
        summary['Gaps'] = MSA_gaps.iloc[:, 0]
        summary['Insertions'] = MSA_insertions.iloc[:, 0]
        summary['Deletions'] = MSA_deletions.iloc[:, 0]

        summary['CDS_Similarity'] = CDS_similarity.iloc[:, 0]
        summary['CDS_Mutations'] = CDS_mutation_num.iloc[:, 0]
        summary['CDS_Mutation_Freq'] = summary['CDS_Mutations']/CDS_length
        summary['CDS_Transitions'] = CDS_transition_num.iloc[:, 0]
        summary['CDS_Transversions'] = CDS_transversion_num.iloc[:, 0]
        summary['CDS_TT_ratio'] = CDS_tt_ratio.iloc[:, 0]
        summary['CDS_Gaps'] = CDS_gaps.iloc[:, 0]
        summary['CDS_Insertions'] = CDS_insertions.iloc[:, 0]
        summary['CDS_Deletions'] = CDS_deletions.iloc[:, 0]

        summary['nonCDS_Similarity'] = nonCDS_similarity.iloc[:, 0]
        summary['nonCDS_Mutations'] = nonCDS_mutation_num.iloc[:, 0]
        summary['nonCDS_Mutation_Freq'] = summary['nonCDS_Mutations']/nonCDS_length
        summary['nonCDS_Transitions'] = nonCDS_transition_num.iloc[:, 0]
        summary['nonCDS_Transversions'] = nonCDS_transversion_num.iloc[:, 0]
        summary['nonCDS_TT_ratio'] = nonCDS_tt_ratio.iloc[:, 0]
        summary['nonCDS_Gaps'] = nonCDS_gaps.iloc[:, 0]
        summary['nonCDS_Insertions'] = nonCDS_insertions.iloc[:, 0]
        summary['nonCDS_Deletions'] = nonCDS_deletions.iloc[:, 0]

        summary_pairwise = pd.DataFrame(columns=['Similarity', 'Mutations', 'Mutation_Freq', 'Transitions', 'Transversions', 'TT_ratio', 'Gaps', 'Insertions', 'Deletions',
                                                 'CDS_Similarity', 'CDS_Mutations', 'CDS_Mutation_Freq', 'CDS_Transitions', 'CDS_Transversions', 'CDS_TT_ratio', 'CDS_Gaps', 'CDS_Insertions', 'CDS_Deletions',
                                                 'nonCDS_Similarity', 'nonCDS_Mutations', 'nonCDS_Mutation_Freq', 'nonCDS_Transitions', 'nonCDS_Transversions', 'nonCDS_TT_ratio', 'nonCDS_Gaps', 'nonCDS_Insertions', 'nonCDS_Deletions'])

        for i in range(len(names)-1):
            summary_pairwise.at[names[i+1],
                                'Similarity'] = MSA_similarity.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'Mutations'] = MSA_mutation_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'Mutation_Freq'] = summary_pairwise.at[names[i+1], 'Mutations']/total_length
            summary_pairwise.at[names[i+1],
                                'Transitions'] = MSA_transition_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'Transversions'] = MSA_transversion_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'TT_ratio'] = MSA_tt_ratio.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'Gaps'] = MSA_gaps.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'Insertions'] = MSA_insertions.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'Deletions'] = MSA_deletions.loc[names[i+1], names[i]]

            summary_pairwise.at[names[i+1],
                                'CDS_Similarity'] = CDS_similarity.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_Mutations'] = CDS_mutation_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_Mutation_Freq'] = summary_pairwise.at[names[i+1], 'CDS_Mutations']/CDS_length
            summary_pairwise.at[names[i+1],
                                'CDS_Transitions'] = CDS_transition_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_Transversions'] = CDS_transversion_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_TT_ratio'] = CDS_tt_ratio.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_Gaps'] = CDS_gaps.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_Insertions'] = CDS_insertions.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'CDS_Deletions'] = CDS_deletions.loc[names[i+1], names[i]]

            summary_pairwise.at[names[i+1],
                                'nonCDS_Similarity'] = nonCDS_similarity.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'nonCDS_Mutations'] = nonCDS_mutation_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1], 'nonCDS_Mutation_Freq'] = summary_pairwise.at[names[i+1],
                                                                                          'nonCDS_Mutations']/nonCDS_length
            summary_pairwise.at[names[i+1],
                                'nonCDS_Transitions'] = nonCDS_transition_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'nonCDS_Transversions'] = nonCDS_transversion_num.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'nonCDS_TT_ratio'] = nonCDS_tt_ratio.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'nonCDS_Gaps'] = nonCDS_gaps.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'nonCDS_Insertions'] = nonCDS_insertions.loc[names[i+1], names[i]]
            summary_pairwise.at[names[i+1],
                                'nonCDS_Deletions'] = nonCDS_deletions.loc[names[i+1], names[i]]

        # """ Output results to files """
        # TODO - add all files to the output
        # name = file_MSA.replace('input/','').split('.')[0]
        # MSA_mutations.to_csv('output/'+name+'_Mutations.csv')
        # MSA_similarity.to_csv('output/'+name+'_Similarity.csv')
        # MSA_mutation_num.to_csv('output/'+name+'_Mutation_number.csv')
        # MSA_tt_ratio.to_csv('output/'+name+'_TT_ratio.csv')
        # MSA_gaps.to_csv('output/'+name+'_Gaps.csv')
        # MSA_insertions.to_csv('output/'+name+'_Insertions.csv')
        # MSA_deletions.to_csv('output/'+name+'_Deletions.csv')

        # CDS_mutations.to_csv('output/'+name+'_CDS_Mutations.csv')

        summary.to_csv('output/'+file_name+'_summary.csv')
        summary_pairwise.to_csv('output/'+file_name+'_summary_pairwise.csv')

        '''
        
        """
        ======================================================
        Protein Stuff
        ======================================================
        """
        
        Gen_code = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K','AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L','CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q','CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V','GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E','GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S','TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_','TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        
        ### TODO
        # implement the checking of codons
        # Two problems:
        # sometimes CDS isn't divisible by 3, so the actual coding sequence is shorter
        # I've partialy solved it by just removing those not divisible by three (since most of them are)
        # Secondly, some sequences are in the opposite direction, so it has to be taken into account
        # when calculating in which position in the codon is the mutation...
        
        CDS_start_copy = CDS_start.copy()
        CDS_stop_copy = CDS_stop.copy()
        CDS_strand_copy = CDS_strand.copy()
        
        for i in range(len(CDS)):
            if (CDS_stop[i] - CDS_start[i]+1) % 3 != 0:
                CDS_start_copy = CDS_start_copy.drop(i)
                CDS_stop_copy = CDS_stop_copy.drop(i)
                CDS_strand_copy = CDS_strand_copy.drop(i)
                
        CDS_array = set()
        
        for i in CDS_start_copy.index:
            CDS_array = CDS_array.union(set(range(CDS_start_copy[i],CDS_stop_copy[i]+1)))
        
        CDS_df = MSA_df[MSA_df.iloc[:,0].isin(CDS_array)]
        
        # make a dataframe which will only contain "mutated" segments of coding sequences
        CDS_mutations = CDS_df[CDS_df["alignment"]==' '].iloc[:,:-1]
        
        CDS_df_lite = CDS_df.iloc[:,range(1,CDS_df.shape[1],2)]
        
        length_CDS = len(CDS_df_lite)
        
        for name1 in names:
            seq1 = CDS_df_lite[name1].values
            for name2 in names:
                seq2 = CDS_df_lite[name2].values
                synonymous = 0
                nonsynonymous = 0
                Codon_pos1 = 0
                Codon_pos2 = 0
                Codon_pos3 = 0
                for i in range(0,length_CDS,3):
                    codon1 = seq1[i:i+3]
                    codon2 = seq2[i:i+3]
        
        '''
