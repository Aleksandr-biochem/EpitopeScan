import os
import sys
import time
import glob
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
    
    stat = df.iloc[:, 3:-4].apply(lambda x: (x == 1).sum())
    
    # create stat dictionary
    AA_mutations_stat = dict()
    for item in stat.items():
        if 'Δ' in item[0]:
            score = -6.0
        else:
            score = blosum_matrix[item[0][0]][item[0][-1]]
        AA_mutations_stat[item[0]] = (item[1], score)

    return AA_mutations_stat
    
def CombinationsStat(df, blosum_version=90):
    """
    Count and score protein combinationed mutations in data
    
    df - pd DataFrame, mutations data

    Returns:
    AA_mutations_stat - dict(str(mutation): (int(count), float(score)))
    """
    blosum_matrix = bl.BLOSUM(blosum_version)

    # create the stat dictionary
    AA_mutations_stat = dict()
    
    # get mutation combinations
    all_mutations = set(df.columns[3:-4])
    mutations, combinations = df.iloc[:, 3:-4], []
    if mutations.shape[1] > 0:
        combinations = mutations.drop_duplicates()
        combinations = combinations.apply(lambda row: list(row[row == 1].index), axis=1)

    # count and score combinations
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
                    score += blosum_matrix[mut[0]][mut[-1]]

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

    Returns:
    summary_table - pd DataFrame
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

    # remove 0 counts
    AA_mutations_stat = {k:v for k, v in AA_mutations_stat.items() if v[0] > 0}
    
    summary_table = []
    for mut in AA_mutations_stat:
        
        summary_row = []
        mut_string = ['_'] * len(peptide)

        for m in mut.split(','):

            position = int(m[1:-1]) 
            AA2 = '-' if ('Δ' in m) else m[-1]
            mut_string[position - peptide.protein_start] = AA2

        summary_row.append(''.join(mut_string))
        summary_row.append(f"{mut}")
        summary_row.append(f"{AA_mutations_stat[mut][0]}")
        summary_row.append(f"{AA_mutations_stat[mut][1]}")
        
        summary_table.append(summary_row)
    
    # print the summary
    summary_table = pd.DataFrame(summary_table, columns = [peptide.sequence, ' ', 'count', 'score'])

    return summary_table

def StatEpitopeData(input_dir,
					proteome,
					reference_genome,
                    stat_mutations=0,
                    sort_by=0,
                    time_start=None,
                    time_end=None,
                    blosum_version=90,
                    metadata_filter=False):
    """
    Print mutations summary from EpitopeScan results data

    input_dir - str, path to EpitopeScan results dir
    proteome - list(Protein instances), reference proteome to map peptides
    reference_genome - str, reference genome sequence
    stat_mitations - int, 0 (stat individual mutations) or 1 (combinations)
    sort_by - int, 0 (by count) or 1 (by BLOSUM score)
    time_start - datetime, start date to stat mutations
    end_time - datetime, end date to stat mutations
    blosum_version - int, BLOSUM version to use for mutation scoring
    metadata_filter - bool, filter samples with metadata
    """

    # check if we work with one or several peptides
    list_dir = [os.path.isfile(os.path.join(input_dir, f)) for f in os.listdir(input_dir)]
    one_peptide = all(list_dir)

	# list peptide names to analyse 
    if one_peptide:
        peptide_names = glob.glob(f"{input_dir}/*.tsv")[0]
        peptide_names = [peptide_names.split('/')[-1].split('_')[0]]
        prefixes = [f"{input_dir}/"]
    else:
        peptide_names = next(os.walk(input_dir))[1]
        prefixes = [f"{input_dir}/{name}/" for name in peptide_names]

    # read data and calculate statistics for each peptide
    printed_headline_statistics, message = False, ''
    for name, prefix in zip(peptide_names, prefixes):

        # extract peptide sequence
        matrix = pd.read_csv(f"{prefix}{name}_AA_mutation_matrix.csv", index_col=0)
        peptide = Protein(name, ''.join(matrix.index))

        # map peptide onto proteome
        peptide = MapPeptides([peptide],
                              proteome,
                              reference_genome,
                              verbose=False)[0]

        # read mutation data
        df = pd.read_csv(f"{prefix}{name}_mutation_data.tsv",
                         sep='\t', index_col=0)

        # convert sample_date to datetime
        df['sample_date']= pd.to_datetime(df['sample_date'], format='%Y-%m-%d')

        # report total number of samples
        if not printed_headline_statistics:
            message += f"Found data on {len(peptide_names)} peptides from {df.shape[0]} samples\n"

        # filter samples with metadata if required
        if metadata_filter:
            df = df[df['has_metadata'] == 1]
            if not printed_headline_statistics:
                message += f"Only keeping {df.shape[0]} samples with metadata\n"
        else:
            # in case dates are provided - ignore them
            start_date, end_date = None, None

        # filter by date and get reported range
        if not time_start is None:
            df = df[df['sample_date'] >= time_start]
            start_date = time_start.strftime('%d/%m/%Y')
        else:
            start_date = df['sample_date'].min().strftime('%d/%m/%Y')

        if not time_end is None:
            df = df[df['sample_date'] <= time_end]
            end_date = time_end.strftime('%d/%m/%Y')
        else:
            end_date = df['sample_date'].max().strftime('%d/%m/%Y')

        # headline statistics
        if not printed_headline_statistics:
            message += f"{df.shape[0]} samples dated between {start_date} and {end_date}\n"
            print(message)
            printed_headline_statistics = True

        # print corresponding summary
        summary = MutationTextSummary(peptide,
		                            df,
		                            stat_mutations,
		                            sort_by,
		                            blosum_version)
        
        print(print(f"{summary.shape[0]} mutations for {peptide.name}"))
        print(summary)
        print()

    return