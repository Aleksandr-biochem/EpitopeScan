import os
import sys
import time
import blosum as bl
import pandas as pd

sys.path.append(os.path.realpath(os.path.dirname(__file__)))
from ProteinUtils import *

def StatIndividualMutations(df, blosum_version=90):
    """
    Count and score individual protein mutations in data
    
    df - pd DataFrame, mutations data

    Returns:
    AA_mutations_stat - dict(str(mutation): (int(count), float(score)))
    """
    blosum_matrix = bl.BLOSUM(blosum_version)
    
    stat = df.iloc[:, 3:].apply(lambda x: (x == 1).sum())
    
    # create stat dictionary
    AA_mutations_stat = dict()
    for item in stat.items():
        if 'Δ' in item[0]:
            score = -6.0
        else:
            score = blosum_matrix[f"{item[0][0]}{item[0][-1]}"]
        AA_mutations_stat[item[0]] = (item[1], score)

    return AA_mutations_stat
    
def CombinationsStat(df, blosum_version=90):
    """
    Count and score protein mutations combinations in data
    
    df - pd DataFrame, mutations data

    Returns:
    AA_mutations_stat - dict(str(mutation): (int(count), float(score)))
    """
    blosum_matrix = bl.BLOSUM(blosum_version)

    # create the stat dictionary
    AA_mutations_stat = dict()
    
    # get mutation combinations
    all_mutations = set(df.columns[3:])
    mutations, combinations = df.iloc[:, 3:], []
    if mutations.shape[1] > 0:
        combinations = mutations.drop_duplicates()
        combinations = combinations.apply(lambda row: list(row[row == 1].index), axis=1)

    # count and scor ecombinations
    for combination in combinations:
        if len(combination) > 0:
            
            mutations_filt = mutations

            for mut in combination:
                mutations_filt = mutations_filt[mutations_filt[mut] == 1]
            for mut in (all_mutations - set(combination)):
                mutations_filt = mutations_filt[mutations_filt[mut] == 0]

            count = mutations_filt.shape[0]

            score = 0.0
            for mut in combination:
                if 'Δ' in mut:
                    score -= 6.0
                else:
                    score += blosum_matrix[f"{mut[0]}{mut[-1]}"]

            AA_mutations_stat[','.join(combination)] = (count, score)
    
    return AA_mutations_stat

def MutationTextSummary(peptide, df,
                        stat_mutations=0,
                        sort_by=0,
                        blosum_version=90):
    """
    Print protein mutation summary to stdout

    peptide - Protein instance, peptide to summarise
    df - pd DataFrame, mutation data
    stat_mutations - int, 0 to stat individual mutations, 1 - mutation combinations
    sort_by - int, 0 to sort mutation report by count, 1 - by BLOSUM score
    blosum_version - BLOSUM version to use for scoring, default 90
    """

    # count and score mutations
    if stat_mutations == 0:
        AA_mutations_stat = StatIndividualMutations(df, blosum_version)
    elif stat_mutations == 1:
        AA_mutations_stat = CombinationsStat(df, blosum_version)
    
    # sort AA substitutions by counts or scores
    key_to_sort = {0: (lambda x: x[1][0]), # by count 
                   1: (lambda x: x[1][1])} # by score
    AA_mutations_stat = dict(sorted(AA_mutations_stat.items(),
                                    key=key_to_sort[sort_by],
                                    reverse=True))
    
    print(f"{len(AA_mutations_stat)} mutations for {peptide.name}")
    print(peptide.sequence)

    for mut in AA_mutations_stat:
        
        summary_string = ['_'] * len(peptide)

        for m in mut.split(','):

            position = int(m[1:-1]) 
            AA2 = '-' if ('Δ' in m) else m[-1]
            summary_string[position - peptide.protein_start] = AA2

        summary_string = ''.join(summary_string)
        summary_string += f"\tcount: {AA_mutations_stat[mut][0]}"
        summary_string += f"\tscore: {AA_mutations_stat[mut][1]}"
        
        print(summary_string)
    
    print()                    
    return

def StatEpitopeData(input_dir, proteome,reference_genome,
                    stat_mutations=0, sort_by=0,
                    # time_start=None, time_end=None,
                    blosum_version=90):
    """
    Print statistics from EpitopeScan results data

    input_dir - str, path to EpitopeScan results dir
    """

    # list peptide names to analyse 
    peptide_names = next(os.walk(input_dir))[1]

    # read data and calculate statistics for each peptide
    printed_headline_statistics = False
    for name in peptide_names:

        # extract peptide sequence
        matrix = pd.read_csv(f"{input_dir}/{name}/{name}_AA_mutation_matrix.csv", index_col=0)
        peptide = Protein(name, ''.join(matrix.index))

        # map peptide onto proteome
        peptide  = MapPeptides([peptide],
                                proteome,
                                reference_genome,
                                verbose=False)[0]

        # read mutation data
        df = pd.read_csv(f"{input_dir}/{name}/{name}_mutation_data.tsv",
                         sep='\t', index_col=0)

        # headline statistics
        if not printed_headline_statistics:
            print(f"Found data on {len(peptide_names)} epitopes from {df.shape[0]} samples\n")
            printed_headline_statistics = True

        # print corresponding summary
        MutationTextSummary(peptide,
                            df,
                            stat_mutations,
                            sort_by,
                            blosum_version)

    return