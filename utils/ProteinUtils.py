class Protein:
    """
    Class Protein

    Attributes:
    name - str, protein name
    sequence - str, protein sequence 
    genome_start - int, start coordinate in reference genome
    genome_end - int, end coordinate in reference genome
    parent_protein - list(str), names of parent proteins
    coding_sequence - str, coding DNA sequence 
    protein_start - int, start coordinate in a primary parent protein
    AA_mutations_matrix - np array, AA substitutions and deletions counts
    NA_mutations_matrix - np array, NA substitutions and deletions counts
    """
    
    def __init__(self,
                 name,
                 sequence):
        
        self.name = name
        self.sequence = sequence

        self.genome_start = None
        self.genome_end = None

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
        Locate peptide among the list of proteins
        Assign genome coordinates accordingly

        proteins - list(Protein instances), proteins to scan
        """

        ###############################################
        # 1st Plp1ab resisue produced after frame shift
        prot_shift_coord = 4402
        ###############################################

        for protein in proteins:
                
            # scan the protein to locate the sequence
            for i in range(0, len(protein)-len(self)):

                # slice subsequence to compare
                sebsequence = protein[i:i + len(self)]

                if sebsequence == self.sequence:

                    # calculate genome coordinates for the peptide
                    genome_start_coordinate = protein.genome_start + (i * 3)
                    genome_end_coordinate = genome_start_coordinate + (len(self) * 3) - 1

                    # account for shift in Plp1ab if any
                    if (protein.name == 'Plp1ab'):
                        if i >= prot_shift_coord - 1:
                            genome_start_coordinate -= 1
                            genome_end_coordinate -= 1
                        else:
                            if (i + len(self) - 1) >= prot_shift_coord - 1:
                                genome_end_coordinate -= 1
                    
                    # assign found coordinates
                    if self.genome_start is None:
                        self.genome_start = genome_start_coordinate
                        self.genome_end = genome_end_coordinate
                        self.protein_start = i + 1
                        self.parent_protein.append(protein.name)

                    # check the case of multiple match
                    else:
                        # in case of same genome coordinate for overlapping ORFs
                        if self.genome_start == genome_start_coordinate:
                            self.parent_protein.append(protein.name)
                        # in case of ambiguous location
                        else:
                            self.genome_start = -1
                            self.genome_end = -1
                            if verbose:
                                print(f"{self.name} has ambiguous location.")
                                print(f"First located at {self.genome_start}({self.parent_protein[-1]}), then at {protein.genome_start + (i * 3)}({protein.name})")
        
        return self

    def AssignCodingSequence(self, reference_genome):
        """
        Assign reference coding DNA sequence
        based on genome start and end coordinates

        reference_genome - str, reference genome sequence
        """

        if (not self.genome_start is None) and (not self.genome_start == -1):
            self.coding_sequence = reference_genome[self.genome_start - 1:self.genome_end]
        else:
            peint(f"Warning! Cannot assign coding sequence to {self.name} due to indefinite coordinates")
        
        return self

def ReadProteinsFromFile(file, proteome=None):
    """
    Read proteins from fasta file

    file - str, path to fasta file
    proteome - list(Protein instances), reference proteome (optional)

    Return:
    list(Protein instances)
    """
    proteins = []

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name_data = line[1:].strip()
                genome_start, genome_end = None, None

                # get start coordinate for reference proteins
                if '/' in name_data:
                    name, genome_start = name_data.split('/')
                    genome_start = int(genome_start)

                # get parent protein and coordinates for special input
                elif ',' in name_data:
                    name_data = name_data.split(',')
                    if len(name_data) != 4:
                        raise Exception(f"expected Name,Patent protein,start, end got: {','.join(name_data)}")
                    else:
                        # get the sequence
                        name, parent_protein, start_res, end_res = name_data
                        start_res, end_res = int(start_res), int(end_res)
                        seq = ''
                        for protein in proteome:
                            if protein.name == parent_protein:
                                seq = protein[start_res - 1:end_res]
                        if seq != '':
                            new_protein = Protein(name, seq)
                        else:
                            raise Exception(f"Couldnot locate peptide {','.join(name_data)}. Check parent protein name and coordinates")
                        proteins.append(new_protein)
                
                # just peptide name followed by sequence
                else:
                    name = name_data

            else:
                seq = line.strip()
                if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
                    print(f"Warning! Protein {name} sequence contains unrecognised characters")
                else:
                    if not genome_start is None:
                        genome_end = genome_start + (len(seq) * 3) - 1
                    new_protein = Protein(name, seq)
                    new_protein.genome_start = genome_start
                    new_protein.genome_end = genome_end
                    proteins.append(new_protein)

    return proteins

def LoadPeptideInput(single_epitope, epitopes_file, proteome, verbose=True):
    """
    Check input and load peptide epitopes
    into Protein instances

    single_epitope - str, single peptide specification
    epitopes_file - str, file with multiple epitopes
    proteome - list(Protein instances), reference proteome

    Returns:
    list(Protein instances) - peptides loaded from input
    """
    # for protein in proteome: 
    #     print(protein.name)

    # read single input peptide
    if not single_epitope is None:
        peptide_data = single_epitope.split(',')

        # name + sequence
        if len(peptide_data) == 2:
            name, seq = peptide_data
            if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
                raise Exception("Peptide sequence contains unrecognised characters")
            epitopes_to_scan = [Protein(name, seq)]
        
        # parent protein + coodinates
        elif len(peptide_data) == 4:
            name, parent_protein, start_res, end_res = peptide_data
            start_res, end_res = int(start_res), int(end_res)

            # get the sequence
            seq = ''
            for protein in proteome:
                if protein.name == parent_protein:
                    seq = protein[start_res - 1:end_res]
            if seq != '':
                epitopes_to_scan = [Protein(name, seq)]
            else:
                raise Exception("Could not locate peptide. Check parent protein name and coordinates")

        else:
            raise Exception("Couldnot parse peptide input. Check format")

    # read multiple peptides from input file
    elif not epitopes_file is None:
        epitopes_to_scan = ReadProteinsFromFile(epitopes_file, proteome)
        if len(epitopes_to_scan) == 0:
            raise Exception("Could not recognise any peptides from input file")

    if verbose:
        print("Input epitopes:")
        for epitope in epitopes_to_scan:
            print(f"{epitope.name} {epitope.sequence}")
        print()

    return epitopes_to_scan

def MapPeptides(peptides, proteome, ref_genome, verbose=True):
    """
    Map peptides onto reference proteome and assign coding sequences

    peptides - list(Protein instances), peptides to map
    proteome - list(Protein instances), reference proteins to map onto
    ref_genome - str, reference genome sequence
    verbose - bool, print mapping statistics, default True

    Returns:
    list(Protein instances), list of unambiguously mapped Protein instances
    """

    couldnt_map, mapped_ambiguously = 0, 0
    mapped_peptides = [] # list to append mapped peptides

    for peptide in peptides:

        peptide.LocateSubsequence(proteome)

        # check if peptide was located successfully
        if peptide.genome_start is None:
            couldnt_map += 1
        # check ambiguous location
        elif peptide.genome_start == -1:
            mapped_ambiguously += 1
        else:
            # assign coding DNA sequence for mapped peptide
            peptide.AssignCodingSequence(ref_genome)
            mapped_peptides.append(peptide)

    if len(mapped_peptides) == 0:
        raise Exception('Could not map any peptides')

    if verbose:
        print(f"Mapped {len(mapped_peptides)} peptides")
        print(f"Could not map {couldnt_map}")
        print(f"Mapped ambiguously {mapped_ambiguously}\n")

    return mapped_peptides






