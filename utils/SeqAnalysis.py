import os
import re
import time
import numpy as np
import blosum as bl
import pandas as pd

def translate(seq, codon_table):
    """
    Translate DNA sequence into amino acid sequence
    indicating deletions and ambiguous residues
    assumes input sequence of 3n

    seq - str, DNA sequence of length 3n
    codon_table - dict, DNA codon table

    Returns:
    str - translated protein sequence
    """

    protein = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon == '---':
            protein += '-'
        elif codon not in codon_table.keys():
            protein += 'X'
        else:
            protein += codon_table[codon]

    return protein

def GlobalPtoteinAlignment(seq1, seq2, blosum_version=90, gap_penalty=-2.0):
    """
    Needleman–Wunsch algorithm alignment for 2 protein sequences 
    with BLOSUM scoring matrix
    Implemetation based on:
    https://python.plainenglish.io/global-sequence-alignment-implementation-in-python-from-scratch-5a8a611dbb1e
    
    seq1 - str, 1st protein sequence
    seq2 - str, 2nd protein sequence
    blosum_version - int, BLOSUM version for scoring
    gap_penalty - float, gap penalty for alignment
    
    Returns:
    (seq1_aligned, seq1_aligned) - (str, str) aligned protein sequences
    """
    BLOSUM_matrix = bl.BLOSUM(blosum_version)
    
    # create intial matrix
    scoring_matrix = np.zeros((len(seq1)+1, len(seq2)+1))
    # initialise first column and first row
    for i in range(1, len(seq1)+1):
        scoring_matrix[i][0] = scoring_matrix[0][0] + (i * gap_penalty)
    for i in range(1, len(seq2)+1):
        scoring_matrix[0][i] = scoring_matrix[0][0] + (i * gap_penalty)
        
    # initialise trace-back matrix
    trace_back_matrix = np.zeros((len(seq1)+1, len(seq2)+1), dtype=str)
    for i in range(1, len(seq1)+1):
        trace_back_matrix[i][0] = 'u'
    for i in range(1, len(seq2)+1):
        trace_back_matrix[0][i] = 'l'
        
    # fill the score and trace-back matrices
    for i in range(1, len(seq1)+1): # rows
        for j in range(1, len(seq2)+1): # columns
            
            # define score for each transition
            left = scoring_matrix[i][j - 1] + gap_penalty
            up   = scoring_matrix[i - 1][j] + gap_penalty

            # for ambiguous bases assume that it is the same one
            if seq2[j-1] == 'X':
                score = BLOSUM_matrix[f"{seq1[i-1]}"][f"{seq1[i-1]}"]
            else:
                score = BLOSUM_matrix[f"{seq1[i-1]}"][f"{seq2[j-1]}"]

            diag = scoring_matrix[i - 1][j - 1] + score
            
            # define the maximum score
            scoring_matrix[i][j] = max(left, up, diag)
            
            # fill traceback
            if scoring_matrix[i][j] == left:
                trace_back_matrix[i][j] = 'l'
            elif scoring_matrix[i][j] == up:
                trace_back_matrix[i][j] = 'u'
            else:
                trace_back_matrix[i][j] = 'd'
    
    # get trace back to produce alignment
    seq1_aligned, seq2_aligned = [], []
    i, j = len(seq1), len(seq2)
    while (i > 0 or j > 0):
        if trace_back_matrix[i][j] == 'd':
            seq1_aligned.append(seq1[i-1])
            seq2_aligned.append(seq2[j-1])
            i -= 1
            j -= 1
        elif trace_back_matrix[i][j] == 'l':
            seq1_aligned.append('-')
            seq2_aligned.append(seq2[j-1])
            j -= 1
        elif trace_back_matrix[i][j] == 'u':
            seq1_aligned.append(seq1[i-1])
            seq2_aligned.append('-')
            i -= 1
        elif trace_back_matrix[i][j] == '':
            break
            
    seq1_aligned = ''.join(seq1_aligned[::-1])
    seq2_aligned = ''.join(seq2_aligned[::-1])
    
    return seq1_aligned, seq2_aligned

def CompareNAsequences(peptide, sample_seq):
    """
    Compare two DNA sequences to report substitutions and deletions
    
    peptide - Protein instance, reference peptide for comparison
    sample_seq - str, sample DNA sequence to compare

    Returns
    list(str) list of mutations
    """
    
    NA_substitutions = [] # list of NA substitutions and deletions

    # compare all bases
    for i, (ref_base, base) in enumerate(zip(peptide.coding_sequence, sample_seq)):

        # report deletion
        if base == '-':
            NA_substitutions.append(f"Δ{peptide.genome_start + i}{ref_base}")

        # report substitution for unambiguous bases
        elif (ref_base != base) and (base in 'ATCG'):
            NA_substitutions.append(f"{ref_base}{peptide.genome_start + i}{base}")

        # count unambiguous bases and deletions
        if base in 'ATCG-':
            ind_NA = 'ATCG-'.find(base)
            peptide.NA_mutations_matrix[i, ind_NA] += 1

    return NA_substitutions

def CompareAAsequences(peptide,
                       sample_seq,
                       codon_table,
                       ambiguity_threshold,
                       partially_spanning_deletions_found,
                       add_seq_left, add_seq_right):
    """
    Compare DAN sequences to report AA substitutions and deletions
    
    peptide - Protein instance, peptide for comparison
    sample_seq - str, sample DNA sequence to compare
    codon_table - dict, DNA codon table
    ambiguity_threshold - float, max proportion of ambiguous peptide residues
    in sample sequence to treat as sufficient coverage
    partially_spanning_deletions_found - bool, True if sequence contains
    partially spanning deletions
    add_seq_left - str, sequence to be added to the left edge of reference coding sequence
    add_seq_right - str, sequence to be added to the right edge of reference coding sequence

    Returns:
    list(str) AA substitutions and deletions
    """
    # list of AA mutations
    AA_substitutions = []

    refrnc_seq = peptide.coding_sequence
    # expand reference coding sequence for frame shift
    # genome coordinate of ORF1ab -1 frameshift is 13469
    if (peptide.genome_start < 13469 - 1) and (13469 < peptide.genome_end):
        refrnc_seq = refrnc_seq[:13469 - peptide.genome_start] + \
                     refrnc_seq[(13469 - peptide.genome_start) - 1:]

    # for partially spanning deletions
    if partially_spanning_deletions_found:

        # expand edges if any
        refrnc_seq = add_seq_left + refrnc_seq + add_seq_right

        # translate DNA sequences without deletions
        sample_protein = translate(sample_seq.replace('-', ''), codon_table)
        refrnc_protein = translate(refrnc_seq, codon_table)

        # call global alignent 
        refrnc_protein, sample_protein = GlobalPtoteinAlignment(refrnc_protein, sample_protein)

        # cut aligned sequences if needed
        sample_protein = sample_protein[len(add_seq_left)//3:]
        refrnc_protein = refrnc_protein[len(add_seq_left)//3:]
        if len(add_seq_right) > 0:
            sample_protein = sample_protein[:-(len(add_seq_right)//3)]
            refrnc_protein = refrnc_protein[:-(len(add_seq_right)//3)]
    # if no spanning deletions are found just translate sequences
    else:
        sample_protein = translate(sample_seq, codon_table)
        refrnc_protein = translate(refrnc_seq, codon_table)

    # check content of ambiguous residues
    ambiguous_residues = sample_protein.count('X')
    if (ambiguous_residues / len(sample_protein)) > ambiguity_threshold:
        AA_substitutions = ['-']
    # process sequences with sufficient coverage
    else:
        # compare residue by residue
        for i, (ref_res, sample_res) in enumerate(zip(refrnc_protein, sample_protein)):
            
            coordinate = peptide.protein_start + i

            # record deletion
            if sample_res == '-':
                AA_substitutions.append(f'Δ{coordinate}{ref_res}')
                # count deletion
                peptide.AA_mutations_matrix[i, -1] += 1
           
           # record substitution for unambiguous residues
            elif sample_res != 'X':
                # for unambiguous codons report residue change if any
                if ref_res != sample_res:
                    coordinate = peptide.protein_start + i
                    AA_substitutions.append(f"{ref_res}{coordinate}{sample_res}")

                # count residue
                ind_AA = 'GALMFWKQESPVICYHRNDT*'.find(sample_res)
                peptide.AA_mutations_matrix[i, ind_AA] += 1
    
    return AA_substitutions

def ReportSequenceMuatations(peptide,
                             genome_seq,
                             reference_genome,
                             codon_table,
                             ambiguity_threshold):
    """
    Report NA and AA mutations in peptide sequence 
    against reference sequence

    peptide - Protein instance, peptide for comparison
    genome_seq - str, aligned sample genome sequence
    reference_genome - str, reference viral genome
    codon_table - dict, DNA codon table
    ambiguity_threshold - float, max proportion of ambiguous peptide residues
    in sample sequence to treat as sufficient coverage

    Returns:
    list(str, str) - NA mutations and AA mutations
    """

    ########################################
    # ORF1ab -1 frameshift genome coordinate
    orf_shift_coord = 13469
    ########################################

    # get the coding sequence for NA comparison
    sample_seq_for_NA_comparison = genome_seq[peptide.genome_start - 1:peptide.genome_end]

    # get the coding sequence for AA comparison
    # in case of frame shift extend for codon-by-codon comparison
    if (peptide.parent_protein[0] == 'Plp1ab') and \
       ((peptide.genome_start < orf_shift_coord) and \
       (orf_shift_coord < peptide.genome_end)):
        sample_seq_for_AA_comparison = genome_seq[peptide.genome_start - 1:orf_shift_coord - 1] + \
                                       genome_seq[orf_shift_coord - 2:peptide.genome_end]
    else:
    	sample_seq_for_AA_comparison = genome_seq[peptide.genome_start - 1:peptide.genome_end] 

    # check for partially spanning deletions
    partially_spanning_deletions_found = False
    for i in range(0, len(sample_seq_for_AA_comparison), 3):
        if sample_seq_for_AA_comparison[i:i + 3].count('-') in [1, 2]:
            partially_spanning_deletions_found = True
            break

    # if partially spanning deletions are found
    # check if these deletions are at sequece edges
    add_seq_left, add_seq_right = '', ''
    if partially_spanning_deletions_found:

        # sequence pretreatment for left edge spanning deletions
        add_left = 0 # n codons to be added on the left
        if sample_seq_for_AA_comparison[0] == '-':

            # deduce how long is this deletion
            del_end = re.search(r"^-+", sample_seq_for_AA_comparison).end()
            del_end_j = del_end % 3

            if del_end_j > 0:
                
                # find bases which will shift to there
                pos = peptide.genome_start - 1
                add_left = 1
                while genome_seq[pos - 3:pos] == '---':
                    add_left += 1
                    pos -= 3

                # edit sequence for AA comparison
                sample_seq_for_AA_comparison = genome_seq[pos - 3:pos + (3 * (add_left - 1))] + \
                                               sample_seq_for_AA_comparison

                # define left flank to ref sequence
                add_seq_left = reference_genome[pos - 3:pos + (3 * (add_left - 1))]

        # sequence pretreatment for right edge deletions
        add_right = 0 # n codons to be added on the right
        if sample_seq_for_AA_comparison[-1] == '-':

            # deduce how long is this deletion
            del_start = re.search(r"-+$", sample_seq_for_AA_comparison).start()
            del_start_j = (len(sample_seq_for_AA_comparison) - del_start) % 3

            if del_start_j > 0:
                
                # find bases which will shift to there
                pos = peptide.genome_end
                add_right = 1
                while genome_seq[pos:pos + 3] == '---':
                    add_right += 1
                    pos += 3

                # edit sequence for AA comparison
                sample_seq_for_AA_comparison = sample_seq_for_AA_comparison + \
                                               genome_seq[pos - (3 * (add_right - 1)):pos + 3]

                # define right flank to ref sequence
                add_seq_right = reference_genome[pos - (3 * (add_right - 1)):pos + 3]

    # compare sequences to report AA mutations
    AA_mutations = CompareAAsequences(peptide,
                                      sample_seq_for_AA_comparison,
                                      codon_table,
                                      ambiguity_threshold,
                                      partially_spanning_deletions_found,
                                      add_seq_left, add_seq_right)

    # if insufficient coverage
    if AA_mutations == ['-']:
        NA_mutations = ['-']
    else:
        # compare DNA sequences 
        NA_mutations = CompareNAsequences(peptide,
                                          sample_seq_for_NA_comparison)

    NA_mutations = ','.join(NA_mutations) if len(NA_mutations) > 0 else None
    AA_mutations = ','.join(AA_mutations) if len(AA_mutations) > 0 else None

    return [NA_mutations, AA_mutations]

def ExpandMutationData(dfs):
    """
    Transform AA mutation column in DataFrame One-Hot Encoding style

    dfs - list(pd DataFrame), list of DataFrames to expand

    Returns:
    list(pd DataFrame) - list of expanded DataFrames
    """

    # for each data frame in list
    for i in range(len(dfs)):

        # get unique protein mutations
        unique_AA_mutations = set()
        for muts in dfs[i]['AA_mutations']:
            if (not muts is None) and (not muts in ['-', 'NF']):
                muts_list = muts.split(',')
                for mut in muts_list:
                    unique_AA_mutations.add(mut)

        # if noo AA mutations found
        if len(unique_AA_mutations) == 0:
        	continue

        # create dataframe of zeros with mutations columns
        AA_substitution_df = pd.DataFrame(0,
                                          index=np.arange(dfs[i].shape[0]),
                                          columns=list(unique_AA_mutations))
        
        # fill dataframe columns
        for j, muts in enumerate(dfs[i]['AA_mutations']):
            if (not muts is None) and (not muts in ['-', 'NF']):
                muts_list = muts.split(',')
                for mut in muts_list:
                    AA_substitution_df.loc[j, mut] = 1

        # concat dataframes
        dfs[i] = pd.concat([dfs[i], AA_substitution_df], axis=1)

    return dfs

def BindMetadata(dfs, metadata_file, sample_tag):
    """
    Bind metadata to mutation DataFrames

    dfs - list(pd DataFrame), list of DataFrames to expand
    metadata_file - str, path to metadata csv file
    sample_tag - str, sample tag to subset

    Returns:
    list(pd DataFrame), list of DataFrames with new metadata columns
    """
	
    # compile regex pattern from tag to filter samples if any
    sample_tag = re.compile(sample_tag) if sample_tag else None

    if not metadata_file is None:
        print(f"Merging output with metadata from {metadata_file}...")

        # read metadata and select needed columns
        metadata = pd.read_csv(metadata_file, low_memory=False)
        metadata = metadata[['sequence_name', 'sample_date',
                             'epi_week', 'usher_lineage']]
        # add column indicating presence of metadata
        metadata['has_metadata'] = 1

        # bind each output df to metadata
        for i in range(len(dfs)):
            dfs[i] = dfs[i].merge(metadata, how='left', on='sequence_name')
            
            # fill nans for samples without metadata
            dfs[i]['has_metadata'] = dfs[i]['has_metadata'].fillna(0)

            # transfrom sample_date to datime format to stat in the end
            dfs[i]['sample_date']= pd.to_datetime(dfs[i]['sample_date'],
            											  format='%Y-%m-%d')

        print(f"Metadata binding finished")
        print(f"Metadata found for {sum(dfs[0]['has_metadata'] == 1)}")
        no_metadata = sum(dfs[0]['has_metadata'] == 0)
        print(f"{no_metadata} samples ({round(no_metadata*100/dfs[0].shape[0], 1)}%) lack metadata\n")

    else:
    # if no metadata given - create dummy columns of Nans
    	for i in range(len(dfs)):
    		dfs[i]['has_metadata'] = 0
    		for col in ['sample_date', 'epi_week', 'usher_lineage']:
    			dfs[i][col] = np.nan	

    return dfs

def CheckFrameDisruption(orf_start, orf_end, genome_sequence):
    """
    Check if the ORF in sample genome is disrupted by deletions
    
    orf_start - int, ORF start coordinate in genome
    orf_end - int, ORF end coordinate in genome
    genome_sequence - str, sample genome sequence

    Returns:
    bool, True if frame is disrupted
    """

    # get ORF sequence
    ORF_sequence = genome_sequence[orf_start - 1:orf_end]

    # check start and stop codons mutations
    if ORF_sequence[:3] != 'ATG': 
        frame_disrupted = True
        return frame_disrupted
    if not ORF_sequence[-3:] in ['TAA', 'TGA', 'TAG']:
        frame_disrupted = True
        return frame_disrupted

    # check deletions disrupting the frame
    # we are looking for deletions with length != 3n
    frame_disrupted = False
    if '-' in ORF_sequence:  
        len_gap = 0
        for base in ORF_sequence:
            if base == '-':
                len_gap += 1
            elif len_gap > 0:
                if len_gap % 3 > 0:
                    frame_disrupted = True
                    break
                else:
                    len_gap = 0

    return frame_disrupted

def ScanMSA(epitopes_to_scan,
            msa_file,
            reference_genome,
            sample_tag=None,
            verbose=True,
            quality_filter=None,
            ambiguity_threshold=0.33):
    """
    Scan genome MSA to report mutations within input peptides

    epitopes_to_scan - list(Protein instances), peptides to analyse
    msa_file - str, path to MSA file
    reference_genome - str, reference viral genome
    sample_tag - str, sample tag to subset
    varbose - bool, print key statistics at the end, default False
    quality_filter - float, max N base proportion to tolerate
    ambiguity_threshold - float, max proportion of ambiguous peptide residues
    in sample sequence to treat as sufficient coverage (default 0.33)

    Returns:
    list(pd DataFrame) - list of DataFrames with mutation data for each peptide
    """

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

    total_samples = 0 # total samples count
    rejected_by_tag = 0 # samples rejected by tag
    discarded_by_quality = 0 # count samples discarded by N quality

    # initiate blank substitution matrices for each epitope
    for epitope in epitopes_to_scan:
        epitope.AA_mutations_matrix = np.zeros((len(epitope.sequence), 22), dtype=int)
        epitope.NA_mutations_matrix = np.zeros((len(epitope.coding_sequence), 5), dtype=int)

    # read Open Reading Frames coordinates
    # for ORF disruption control
    cur_dir = os.path.realpath(os.path.dirname(__file__))
    ORF_coordinates = dict()
    with open(f'{cur_dir}/../reference_sequences/ORFs_reference.txt', 'r') as f:
        for line in f:
            orf, start, end, prot, l = line.strip().split()
            ORF_coordinates[prot] = (int(start), int(end))

    # compile regex pattern from input to filter samples if any
    sample_tag = re.compile(sample_tag) if sample_tag else None
    # define if samples have to be filtered by tag
    filter_samples = True if (not sample_tag is None) else False
    
    ## parse MSA file and analyse each sequence
    df = [] # list to append extracted mutatiuon data
    with open(msa_file, 'r') as f:

        # initiate checkpoint time for reports
        time_checkpoint = time.time()

        for line in f:
            # sample name line
            if line.startswith('>'):

                total_samples += 1 # count total parsed samples
                sample_name = line.strip()[1:] # get sample name
                new_row = [sample_name] # initiate new row in mutation data

                # check presense of tag if necessary
                skip_sample = False
                if filter_samples:
                    skip_sample = False if sample_tag.search(sample_name) else True
                    
            # genome sequence line
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

                    # check frame disruptions for epitope parent proteins
                    frame_disruption = dict()
                    for epitope in epitopes_to_scan:
                        if not epitope.parent_protein[0] in frame_disruption.keys():
                            disruption = CheckFrameDisruption(ORF_coordinates[epitope.parent_protein[0]][0],
                                                              ORF_coordinates[epitope.parent_protein[0]][1],
                                                              genome_aln)
                            frame_disruption[epitope.parent_protein[0]] = disruption

                    # compare reference sequence of each epitope to the sample
                    for epitope in epitopes_to_scan:

                        # check if parent protein frame is desrupted 
                        if frame_disruption[epitope.parent_protein[0]]:
                            new_row.append(['NF', 'NF']) # report as non-functional

                        # for non-disrupted translation
                        else:
                            # report NA and AA mutations against reference sequence
                            new_row.append(ReportSequenceMuatations(epitope,
                                                                    genome_aln,
                                                                    reference_genome,
                                                                    codon_table,
                                                                    ambiguity_threshold))

                    # add new row to the data
                    df.append(new_row)
                
                else:
                    # count sample rejected by tag
                    rejected_by_tag += 1

            # update sample count in stdout
            if time.time() - time_checkpoint > 1:
                time_checkpoint = time.time()
                print(f"Scanned samples: {total_samples}", end="\r")
    
    # print scanning statistics
    if verbose:
        print()
        print(f"Found {total_samples} samples")
        if filter_samples:
            print(f"Rejected {rejected_by_tag} samples by tag")
        if not quality_filter is None:
            print(f"Rejected {discarded_by_quality} genomes having >{quality_filter} of N bases")
        print(f"Returning mutation data for {len(df)} samples\n")

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