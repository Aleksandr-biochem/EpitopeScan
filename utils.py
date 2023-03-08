import os
import re
import numpy as np
import blosum as bl
import pandas as pd

class Protein:
    """
    Class Protein

    name - str, protein name
    sequence - str, protein sequence 
    genome_start - int, start coordinate in reference genome
    genome_end - int, end coordinate in reference genome
    parent_protein - list(str), names of parent proteins
    coding_sequence - string, coding DNA sequence 
    protein_start - int, start coordinate in a parent protein
    AA_mutations_matrix - np array, counts of sequence AA substitutions
    NA_mutations_matrix - np array, counts of sequence NA substitutions
    """
    
    def __init__(self,
                 name,
                 sequence,
                 genome_start = None,
                 genome_end   = None):
        
        self.name = name
        self.sequence = sequence
        self.genome_start = genome_start
        if (not genome_start is None) and (genome_end is None):
            # infere genome end coordinate from start and length
            self.genome_end = genome_start + (len(sequence) * 3) - 1
        else:
            self.genome_end = genome_end

        self.parent_protein = list()
        self.protein_start = None
        self.coding_sequence = None
        self.AA_mutations_matrix = None
        self.NA_mutations_matrix = None
        
    def __len__(self):
        return len(self.sequence)
    
    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            item = self.sequence[subscript.start:subscript.stop:subscript.step]
        else:
            item = self.sequence[subscript]
        return item
    
    def __repr__(self):
        return f"{self.name}/{self.sequence}/{len(self)}/{self.genome_start}/{self.genome_end}/{','.join(self.parent_protein)}"
    
    def LocateSubsequence(self, proteins, verbose=False):
        """
        Locate the peptide among list of proteins
        Assign genome coordinate accordingly

        proteins - list(Protein instances), proteins to scan
        """
        for protein in proteins:
                
            # scan the protein to locate the sequence
            for i in range(0, len(protein)-len(self)):

                # slice subsequence to compare
                sebsequence = protein[i:i + len(self)]

                if sebsequence == self.sequence:
                    
                    if self.genome_start is None:
                        self.genome_start = protein.genome_start + (i * 3)
                        self.genome_end = self.genome_start + (len(self) * 3) - 1
                        self.protein_start = i + 1
                        self.parent_protein.append(protein.name)
                    else:
                        # in case of same genome coordinate for overlapping ORFs
                        if self.genome_start == (protein.genome_start + (i * 3)):
                            self.parent_protein.append(protein.name)
                        else:
                            # in case of ambiguous location
                            self.genome_start = -1
                            if verbose:
                                print(f"{self.name} has ambiguous location.")
                                print(f"First located at {self.genome_start}({self.parent_protein[-1]}), then at {protein.genome_start + (i * 3)}({protein.name})")
        
        return

    def AssignCodingSequence(self, reference_genome):
        """
        Assign coding DNA sequence based on start and end coordinates

        reference_genome - str, reference genome sequence
        """
        if (not self.genome_start is None) and (not self.genome_start == -1):
            self.coding_sequence = reference_genome[self.genome_start - 1:self.genome_end]
        return

def ReadProteinsFromFile(file):
    """
    Read protein sequences from fasta file

    file - str, path to fasta file

    Return:
    list(Protein instances)
    """
    proteins = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].strip()
                genome_start = None
                if '/' in line:
                    name, genome_start = name.split('/')
                    genome_start = int(genome_start)
            else:
                seq = line.strip()
                if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
                    print(f"Warning! Protein {name} sequence contains unrecognised characters")
                else:
                    proteins.append(Protein(name, seq, genome_start))
    return proteins

def MapPeptides(peptides, proteome, ref_genome, verbose=True):
    """
    Map peptides onto the proteome and assign coding sequences

    peptides - list of Protein instances, peptides to map
    proteome - list of Protein instances, reference protein to map onto
    ref_genome - str, reference genome sequence
    verbose - bool, print mapping statistics, default True

    Returns:
    mapped_epitopes, list of unambiguously mapped Protein instances
    """
    couldnt_map, mapped_ambiguously = 0, 0
    mapped_peptides = [] # list to append mapped peptides

    for peptide in peptides:

        peptide.LocateSubsequence(proteome)

        # check if peptide was located successfully
        if peptide.genome_start is None:
            couldnt_map += 1
        elif peptide.genome_start == -1:
            mapped_ambiguously += 1
        else:
            # get coding DNA sequence for  successfully mapped peptide
            peptide.AssignCodingSequence(ref_genome)
            mapped_peptides.append(peptide)

    if len(mapped_peptides) == 0:
        raise Exception('Could not map any peptides')

    if verbose:
        print(f"Mapped {len(mapped_peptides)} peptides")
        print(f"Could not map {couldnt_map}")
        print(f"Mapped ambiguously {mapped_ambiguously}\n")

    return mapped_peptides

def CompareSequences(peptide, sample_seq, codon_table):
    """
    Compare DNA sequences to find NA and AA substitutions
    
    peptide - Protein instance, peptide for comparison
    sample_seq - str, sample DNA sequence to compare
    codon_table - dict, DNA codon table

    Returns:
    list( str NA substitutions, str AA substitutions)
    """
    NA_substitutions, AA_substitutions = [], []
    
    # scan DNA codons
    for i in range(0, len(sample_seq), 3):

        ref_codon, codon = peptide.coding_sequence[i:i+3], sample_seq[i:i+3]
        
        # compare bases within the codon
        for j, (ref_base, base) in enumerate(zip(ref_codon, codon)):
            if ref_base != base:
                coordinate = peptide.genome_start + i + j
                NA_substitutions.append(f"{ref_base}{coordinate}{base}")

                # count NA substitution
                ind_NA = 'ATCG'.find(base)
                peptide.NA_mutations_matrix[i + j, ind_NA] += 1
        
        # report residue change if any
        if codon_table[ref_codon] != codon_table[codon]:

            coordinate = peptide.protein_start + (i // 3)
            AA_substitutions.append(f"{codon_table[ref_codon]}{coordinate}{codon_table[codon]}")

            # count AA sunstitution
            ind_AA = 'GALMFWKQESPVICYHRNDT*'.find(codon_table[codon])
            peptide.AA_mutations_matrix[(i // 3), ind_AA] += 1

    NA_substitutions = ','.join(NA_substitutions) if len(NA_substitutions) > 0 else None
    AA_substitutions = ','.join(AA_substitutions) if len(AA_substitutions) > 0 else None    
    
    return [NA_substitutions, AA_substitutions]

def ExpandMutationData(dfs):
    """
    Transform AA mutation column in DataFrame 
    into several columns signifying presence of distinc mutations

    dfs - list(pd DataFrame), list of dataFrames to expand

    Returns:
    list(pd DataFrame) - list of expanded DataFrames
    """

    for i in range(len(dfs)):

        # get unique protein mutations
        unique_AA_mutations = set()
        for muts in dfs[i]['AA_mutations']:
            if (not muts is None) and (muts != '-'):
                muts_list = muts.split(',')
                for mut in muts_list:
                    unique_AA_mutations.add(mut)

        # create dataframe columns
        AA_substitution_df = pd.DataFrame(0,
                                          index=np.arange(dfs[i].shape[0]),
                                          columns=unique_AA_mutations)
        
        # fill dataframe columns
        for j, muts in enumerate(dfs[i]['AA_mutations']):
            if (not muts is None) and (muts != '-'):
                muts_list = muts.split(',')
                for mut in muts_list:
                    AA_substitution_df.loc[j, mut] = 1

        dfs[i] = pd.concat([dfs[i], AA_substitution_df], axis=1)

    return dfs

def ScanMSA(epitopes_to_scan, msa_file, sample_tag=None, verbose=True):
    """
    Scan genome MSA to report mutations within peptides

    epitopes_to_scan - list(Protein instances), peptides to analyse
    msa_file - string, path to MSA file
    sample_tag - compiled regex, sample tag to subset
    varbose - bool, print key statistics at the end, default False

    Returns:
    list(pd DataFrame) - list of DataFrames with mutation data for each peptide
    """

    total_samples = 0 # found samples count

    # initiate blank substitution matrices for each epitope
    for epitope in epitopes_to_scan:
        epitope.AA_mutations_matrix = np.zeros((len(epitope.sequence), 21), dtype=int)
        epitope.NA_mutations_matrix = np.zeros((len(epitope.coding_sequence), 4), dtype=int)
    
    # codon table
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }

    # define if we have to filter samples
    filter_samples = True if (not sample_tag is None) else False
    
    # parse MSA file and analyse each sequence
    df = [] # list to append extracted mutatiuon data
    with open(msa_file, 'r') as f:

        for line in f:

            if line.startswith('>'):

                total_samples += 1
                sample_name = line.strip()[1:]
                new_row = [sample_name]

                # check presense of filtering tag is any
                if filter_samples:
                    skip_sample = False if sample_tag.search(sample_name) else True
                else:
                    skip_sample = False

            else:
                
                if not skip_sample:
                    # get the aligned sequence
                    genome_aln = line.strip()
                    
                    # compare sequence of each epitope to the sample
                    for epitope in epitopes_to_scan:
                        
                        # get DNA sequence corresponding to the epitope
                        sample_sequence = genome_aln[epitope.genome_start - 1:epitope.genome_end]
                        
                        # check that DNA sequence is complete
                        if set(sample_sequence) == set('ATCG'):

                            # run comparison to reference and save mutation data
                            new_row.append(CompareSequences(epitope,
                                                            sample_sequence,
                                                            codon_table))
                        else:
                            # else report missing coverage
                            new_row.append(['-', '-'])

                    # add new row to the data
                    df.append(new_row)
    
    if verbose:
        print(f"Found {total_samples} samples")
        print(f"Rejected {total_samples - len(df)} samples by tag")
        print(f"Returning data for {len(df)} samples\n")

    # separate dataframes for each epitope
    print("Finalising output data...\n")
    dfs = []
    for i in range(1, len(df[0])):
        epitope_data = [[el[0], el[i][0], el[i][1]] for el in df]
        new_df = pd.DataFrame(epitope_data,
                              columns=['sequence_name', 'NA_mutations', 'AA_mutations'])
        dfs.append(new_df)

    # expand DataFrame columns for protein mutations
    dfs = ExpandMutationData(dfs)

    return dfs

def StatIndividualMutations(df, blosum_version=90):
    """
    Count and score individual protein mutations in data
    
    df - pd DataFrame, mutations data

    Returns:
    AA_mutations_stat - dict(str(mutation): (int(count), float(score)))
    """
    blosum_matrix = bl.BLOSUM(blosum_version)
    
    # this is a weak point, explicit specification of columns
    stat = df.iloc[:, 3:].apply(lambda x: (x == 1).sum())
    
    # create stat dictionary
    AA_mutations_stat = dict()
    for item in stat.items():
        AA_mutations_stat[item[0]] = (item[1],
                                      blosum_matrix[f"{item[0][0]}{item[0][-1]}"])
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

            count = mutations_filt.shape[0]
            score = sum([blosum_matrix[f"{mut[0]}{mut[-1]}"] for mut in combination])
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
        
        summary_string = ['-'] * len(peptide)

        for m in mut.split(','):

            position, AA2 = int(m[1:-1]), m[-1]
            summary_string[position - peptide.protein_start] = AA2

        summary_string = ''.join(summary_string)
        summary_string += f"\tcount: {AA_mutations_stat[mut][0]}"
        summary_string += f"\tscore: {AA_mutations_stat[mut][1]}"
        
        print(summary_string)
    
    print()                    
    return

def StatEpitopeData(input_dir,
                    proteome,
                    reference_genome,
                    stat_mutations=0,
                    sort_by=0,
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

        # expract peptide sequence
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





