import os
import re
import time
import numpy as np
import blosum as bl
import pandas as pd


def CompareNAsequences(peptide, sample_seq):
    """
    Compare two DNA sequences to report substitutions and deletions
    
    peptide - Protein instance, peptide for comparison
    sample_seq - str, sample DNA sequence to compare

    Returns
    list(str) list of mutations
    """
    
    NA_substitutions = [] # list of NA sunstitutions and deletions

    # compare all bases
    for i, (ref_base, base) in enumerate(zip(peptide.coding_sequence, sample_seq)):

        # report deletion
        if base == '-':
            NA_substitutions.append(f"Δ{peptide.genome_start + i}{ref_base}")

        # report substitution for unambiguous bases
        elif (ref_base != base) and (base in 'ATCG'):
            NA_substitutions.append(f"{ref_base}{peptide.genome_start + i}{base}")

        # count base change
        if base in 'ATCG-':
            ind_NA = 'ATCG-'.find(base)
            peptide.NA_mutations_matrix[i, ind_NA] += 1

    return NA_substitutions

def CompareAAsequences(peptide, sample_seq, codon_table,
                       ambiguity_intolerance=False):
    """
    Compare DNA sequences to find NA and AA substitutions
    
    peptide - Protein instance, peptide for comparison
    sample_seq - str, sample DNA sequence to compare
    codon_table - dict, DNA codon table
    ambiguity_intolerance - bool, discard any sequence with ambiguous bases

    Returns:
    list(str) AA substitutions and deletions
    """
    AA_substitutions = []
    ambiguous_codons = 0 # report ambiguous codons

    # scan DNA codons
    for i in range(0, len(sample_seq), 3):

        ref_codon, codon = peptide.coding_sequence[i:i+3], sample_seq[i:i+3]

        # deleted codon
        if codon == '---':
            coordinate = peptide.protein_start + (i // 3)
            AA_substitutions.append(f'Δ{coordinate}{codon_table[ref_codon]}')
            peptide.AA_mutations_matrix[i // 3, -1] += 1

        else:

            # report if ambiguous codon is encountered
            if codon not in codon_table.keys():
                ambiguous_codons += 1
            else:

                # for unambiguous codons report residue change if any
                if codon_table[ref_codon] != codon_table[codon]:
                    coordinate = peptide.protein_start + (i // 3)
                    AA_substitutions.append(f"{codon_table[ref_codon]}{coordinate}{codon_table[codon]}")

                # count AA change
                ind_AA = 'GALMFWKQESPVICYHRNDT*'.find(codon_table[codon])
                peptide.AA_mutations_matrix[i // 3, ind_AA] += 1

    # in case of ambiguity intolerance
    if ambiguity_intolerance and ambiguous_codons > 0:
        AA_substitutions = ['-']
    # if too many amb. codons - do not trust
    elif ambiguous_codons > (len(peptide.sequence) // 3):
        AA_substitutions = ['-'] 
    
    return AA_substitutions

def ReportSequenceMuatations(peptide, genome_seq,
                             codon_table, ambiguity_intolerance=False):
    """
    Report NA and AA mutations in peptide sequence 
    against reference sequence

    peptide - Protein instance, peptide for comparison
    genome_seq - str, full aligned sample genome sequence
    codon_table - dict, DNA codon table
    ambiguity_intolerance - bool, discard sequence with ambiguous bases

    Returns:
    list(str, str) - NA mutations and AA mutations
    """

    # get the coding sequence
    sample_seq = genome_seq[peptide.genome_start - 1:peptide.genome_end] 

    # compare DNA and translated protein sequences 

    # if sequence is deleted completely - shortcut 
    if sample_seq.count('-') == len(sample_seq):
        AA_mutations = [f"Δ{peptide.protein_start + i}{r}" for i, r in enumerate(peptide.sequence)]
        peptide.AA_mutations_matrix[:, -1] += 1

    else:

        # sequence pretreatment for left edge deletions
        if sample_seq[0] == '-':

            # deduce how long is this deletion
            del_end = re.search(r"^-+", sample_seq).end()
            del_end_j = del_end % 3

            if del_end_j > 0:
                
                # find bases which will shift to there
                pos = peptide.genome_start - 1
                while genome_seq[pos] == '-':
                    pos -= 1

                # get shifted bases
                shifted_bases = genome_seq[pos - del_end_j + 1:pos + 1]

                # edit sample sequence
                sample_seq = sample_seq[:del_end - del_end_j] + \
                             shifted_bases + \
                             sample_seq[del_end:]

        # sequence pretreatment for right edge deletions
        if sample_seq[-1] == '-':

            # deduce how long is this deletion
            del_start = re.search(r"-+$", sample_seq).start()
            del_start_j = (len(sample_seq) - del_start) % 3

            if del_start_j > 0:
                
                # find bases which will shift to there
                pos = peptide.genome_end - 1
                while genome_seq[pos] == '-':
                    pos += 1

                # get shifted bases
                shifted_bases = genome_seq[pos:pos + del_start_j]

                # edit sample sequence
                sample_seq = sample_seq[:del_start] + \
                             shifted_bases + \
                             sample_seq[del_start + del_start_j:]

        # sequence pretreatment for internal deletions
        gaps = [match.span() for match in re.finditer(r"[A-Z]-+[A-Z]", sample_seq)]
        for gap_start, gap_end in gaps:
            n = 3 - ((gap_start % 3) + 1) # how many bases to move
            sample_seq = sample_seq[:gap_start + 1] + \
                         sample_seq[gap_end - 1:gap_end - 1 + n] + \
                         sample_seq[gap_start + 1:gap_end - 1] + \
                         sample_seq[gap_end - 1 + n:]

        # compare sequences to report mutations
        AA_mutations = CompareAAsequences(peptide,
                                          sample_seq,
                                          codon_table,
                                          ambiguity_intolerance)

    if AA_mutations == ['-']:
        NA_mutations = ['-']
    else:
        # compare DNA sequences 
        NA_mutations = CompareNAsequences(peptide,
                                          genome_seq[peptide.genome_start - 1:peptide.genome_end])

    NA_mutations = ','.join(NA_mutations) if len(NA_mutations) > 0 else None
    AA_mutations = ','.join(AA_mutations) if len(AA_mutations) > 0 else None

    return [NA_mutations, AA_mutations]

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
            if (not muts is None) and (not muts in ['-', 'NF']):
                muts_list = muts.split(',')
                for mut in muts_list:
                    unique_AA_mutations.add(mut)

        # create dataframe columns
        AA_substitution_df = pd.DataFrame(0,
                                          index=np.arange(dfs[i].shape[0]),
                                          columns=unique_AA_mutations)
        
        # fill dataframe columns
        for j, muts in enumerate(dfs[i]['AA_mutations']):
            if (not muts is None) and (not muts in ['-', 'NF']):
                muts_list = muts.split(',')
                for mut in muts_list:
                    AA_substitution_df.loc[j, mut] = 1

        dfs[i] = pd.concat([dfs[i], AA_substitution_df], axis=1)

    return dfs

def CheckFrameDisruption(orf_start, orf_end, genome_sequence):
    """
    Check if the reading frame in sample sequence is disrupted
    
    orf_start - int, ORF start coordinate in genome
    orf_end - int, ORF end coordinate in genome
    genome_sequence - str, sample genome sequence

    Returns:
    frame_disrupted - bool, True is frame is disrupted
    """

    # get ORF sequence
    ORF_sequence = genome_sequence[orf_start - 1:orf_end]

    if '-' in ORF_sequence:  
        frame_disrupted = False
        gap_start, gap_end = None, None

        # scan sequence by codons
        for i in range(0, len(ORF_sequence), 3):
            
            codon = ORF_sequence[i:i+3]
            
            for j, b in enumerate(codon):
                if b == '-':
                    if gap_start is None:
                        gap_start = j
                    gap_end = j
                else:
                    if not gap_start is None:
                        if not (gap_start, gap_end) in [(0, 2), (1, 0), (2, 1)]:
                            frame_disrupted = True
                            return frame_disrupted
                        gap_start, gap_end = None, None
    else:
        frame_disrupted = False

    return frame_disrupted

def ScanMSA(epitopes_to_scan, msa_file, sample_tag=None,
            verbose=True, quality_filter=None,
            ambiguity_intolerance=False):
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
    rejected_by_tag = 0 # samples rejected by tag
    discarded_by_quality = 0 # count samples discarded by quality if neccessary

    # initiate blank substitution matrices for each epitope
    for epitope in epitopes_to_scan:
        epitope.AA_mutations_matrix = np.zeros((len(epitope.sequence), 22), dtype=int)
        epitope.NA_mutations_matrix = np.zeros((len(epitope.coding_sequence), 5), dtype=int)

    # read ORF coordinates for ORF control
    cur_dir = os.path.realpath(os.path.dirname(__file__))
    ORF_coordinates = dict()
    with open(f'{cur_dir}/../reference_sequences/ORFs_reference.txt', 'r') as f:
        for line in f:
            orf, start, end, prot, l = line.strip().split()
            ORF_coordinates[prot] = (int(start), int(end))
    
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

    # define samples have to be filtered by tag
    filter_samples = True if (not sample_tag is None) else False
    
    # parse MSA file and analyse each sequence
    df = [] # list to append extracted mutatiuon data
    with open(msa_file, 'r') as f:

        # record checkpoint time for reports
        time_checkpoint = time.time()

        for line in f:

            if line.startswith('>'):

                total_samples += 1 # count total parssed samples
                sample_name = line.strip()[1:] # get sample name
                new_row = [sample_name]

                # check presense of filtering tag if necessary 
                if filter_samples:
                    skip_sample = False if sample_tag.search(sample_name) else True
                else:
                    skip_sample = False

            else:
                
                if not skip_sample:

                    # get the aligned sequence
                    genome_aln = line.strip()

                    # check proportion of N bases if necessary
                    if not quality_filter is None:
                        proportion_of_N = genome_aln.count('N') / (len(genome_aln) - genome_aln.count('-'))
                        if proportion_of_N > quality_filter:
                            discarded_by_quality += 1
                            continue

                    # check frame disruptions for epitopes parent proteins
                    frame_disruption = dict()
                    for epitope in epitopes_to_scan:
                        if not epitope.parent_protein[0] in frame_disruption.keys():
                            disruption = CheckFrameDisruption(ORF_coordinates[epitope.parent_protein[0]][0],
                                                              ORF_coordinates[epitope.parent_protein[0]][1],
                                                              genome_aln)
                            frame_disruption[epitope.parent_protein[0]] = disruption

                    # compare sequence of each epitope to the sample
                    for epitope in epitopes_to_scan:

                        # check if parent protein frame is desrupted 
                        if frame_disruption[epitope.parent_protein[0]]:
                            new_row.append(['NF', 'NF']) # report as non-functional

                        # for non-disrupted translation
                        else:

                            # report epitope mutations against reference sequence
                            new_row.append(ReportSequenceMuatations(epitope,
                                                                    genome_aln,
                                                                    codon_table,
                                                                    ambiguity_intolerance))

                    # add new row to the data
                    df.append(new_row)
                
                else:
                    rejected_by_tag += 1

            if time.time() - time_checkpoint > 1:
                time_checkpoint = time.time()
                print(f"Scanned samples: {total_samples}", end="\r")
    
    if verbose:
        print()
        print(f"Found {total_samples} samples")
        if filter_samples:
            print(f"Rejected {rejected_by_tag} samples by tag")
        if not quality_filter is None:
            print(f"Rejected {discarded_by_quality} genomes having >0.05 N bases")
        print(f"Returning mutations data for {len(df)} samples\n")

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