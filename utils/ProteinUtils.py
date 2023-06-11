class Protein:
    """
    Class Protein

    Primary attributes:
    name - str, protein name
    sequence - str, protein sequence
    
    Secondary attributes:
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
        
    # allows len(protein)
    def __len__(self):
        return len(self.sequence)
    
    # allows slicing
    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            item = self.sequence[subscript.start:subscript.stop:subscript.step]
        else:
            item = self.sequence[subscript]
        return item

    # defines output when print(protein)
    def __repr__(self):
        return f"{self.name}/{self.sequence}/{len(self)}/{self.genome_start}/{self.genome_end}/{','.join(self.parent_protein)}"
    
    def LocateSubsequence(self, proteins, verbose=False):
        """
        Locate peptide among the list of proteins
        ans assign genome coordinates accordingly

        proteins - list(Protein instances), proteins to scan
        """

        ###############################################
        # 1st Plp1ab resisue produced after frame shift
        prot_shift_coord = 4402
        ###############################################

        for protein in proteins:
                
            # scan the protein to locate the sequence
            for i in range(0, len(protein) - len(self) + 1):

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

    def LocateFromCoordinates(self, parent_protein, start_res, proteome):
        """
        Locate peptide in genome 
        using predefined position in parent protein

        parent_protein - str, predefined parent protein name
        start_res - int, first residue of the peptide in parent protein
        proteome - list(Protein instances), reference proteome
        """

        ###############################################
        # 1st Plp1ab/NSP12 resisue produced after frame shift
        prot_shift_coord = {'Plp1ab': 4402, 'NSP12': 10}
        ###############################################

        # calculate genome coordinates from parent protein
        for protein in proteome:
            if protein.name == parent_protein:
                
                genome_start_coordinate = protein.genome_start + ((start_res - 1) * 3)
                genome_end_coordinate = genome_start_coordinate + (len(self) * 3) - 1
                
                # account for parent proteins containing shift
                if protein.name in prot_shift_coord.keys():
                    if start_res >= prot_shift_coord[protein.name]:
                        genome_start_coordinate -= 1
                        genome_end_coordinate -= 1
                    elif (start_res + len(self) - 1) >= prot_shift_coord[protein.name]:
                            genome_end_coordinate -= 1

        self.genome_start = genome_start_coordinate
        self.genome_end = genome_end_coordinate

        # run through the proteome to define all parent proteins
        # and assign primary protein start
        for protein in proteome:
            if (protein.genome_start <= genome_start_coordinate) and \
               (protein.genome_end   >= genome_end_coordinate):
                
                # for the first primary parent protein
                # calculate protein start minding potential shift
                if len(self.parent_protein) == 0:
                    base_distance = genome_start_coordinate - protein.genome_start
                    if base_distance % 3 == 2:
                        base_distance += 1
                    start = (base_distance // 3) + 1 
                    self.protein_start = start

                # add matching protein to list
                self.parent_protein.append(protein.name)

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
            print(f"Warning! Cannot assign coding sequence to {self.name} due to indefinite coordinates")
        
        return self

def ReadProteinsFromFile(file, proteome=None, ref_genome=None):
    """
    Read proteins from file

    file - str, path to file woth proteins
    proteome - list(Protein instances), reference proteome (optional)
    ref_genome - str, reference genome sequence (optional)

    Return:
    list(Protein instances)
    """
    proteins = [] # list to append proteins

    with open(file, 'r') as f:
        for line in f:

            # FASTA name line or peptide coordinate line
            if line.startswith('>'):
                name_data = line[1:].strip()
                genome_start, genome_end = None, None

                # get start coordinate for reference proteins
                if '/' in name_data:
                    name, genome_start = name_data.split('/')
                    genome_start = int(genome_start)

                # create Protein instance from coordinate input
                elif ',' in name_data:
                    name_data = name_data.split(',')
                    if len(name_data) != 4:
                        raise Exception(f"Expected 'name,parent_protein,start,end'. Got {','.join(name_data)}")
                    else:
                        # get the sequence by coordinates
                        name, parent_protein, start_res, end_res = name_data
                        start_res, end_res = int(start_res), int(end_res)
                        seq = ''
                        for protein in proteome:
                            if protein.name == parent_protein:
                                seq = protein[start_res - 1:end_res]
                        if seq != '':
                            new_protein = Protein(name, seq)
                            new_protein.LocateFromCoordinates(parent_protein,
                                                              start_res,
                                                              proteome)
                        else:
                            raise Exception(f"Couldnot locate peptide {','.join(name_data)}. Check parent protein name and coordinates")
                        
                        proteins.append(new_protein)
                
                # just peptide name followed by sequence
                else:
                    name = name_data

            # sequence line in FASTA file
            else:
                seq = line.strip()
                if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
                    print(f"Warning! Protein {name} sequence contains unrecognised characters")
                else:
                    new_protein = Protein(name, seq)

                    # when reading reference proteins add genome start and end
                    if not genome_start is None:
                        genome_end = genome_start + (len(seq) * 3) - 1
                        new_protein.genome_start = genome_start
                        new_protein.genome_end = genome_end

                    proteins.append(new_protein)

    return proteins

def LoadPeptideInput(single_epitope, epitopes_file, proteome, ref_genome, verbose=True):
    """
    Check input and load peptide epitopes
    into Protein instances

    single_epitope - str, single peptide input
    epitopes_file - str, input file name with multiple epitopes
    proteome - list(Protein instances), reference proteome
    ref_genome - str, reference genome sequence
    verbose - bool, print loaded peptides info (default True)

    Returns:
    list(Protein instances) - peptides loaded from input
    """

    # read single input peptide
    if not single_epitope is None:

        peptide_data = single_epitope.split(',')

        # name + sequence
        if len(peptide_data) == 2:
            name, seq = peptide_data
            # check protein identity of the sequence
            if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
                raise Exception("Peptide sequence contains unrecognised characters")
            epitopes_to_scan = [Protein(name, seq)]
        
        # name + parent protein + coodinates
        elif len(peptide_data) == 4:
            name, parent_protein, start_res, end_res = peptide_data
            start_res, end_res = int(start_res), int(end_res)

            # get the sequence
            seq = ''
            for protein in proteome:
                if protein.name == parent_protein:
                    seq = protein[start_res - 1:end_res]

            # initiate Protein instance with a predefined start
            if seq != '':
                epitopes_to_scan = [
                Protein(name, seq).LocateFromCoordinates(parent_protein,
                                                         start_res,
                                                         proteome)
                ]
            else:
                raise Exception("Could not locate peptide. Check parent protein name and coordinates")

        else:
            raise Exception("Couldnot parse peptide input. Check format")

    # or read multiple peptides from input file
    elif not epitopes_file is None:
        epitopes_to_scan = ReadProteinsFromFile(epitopes_file, proteome, ref_genome)
        if len(epitopes_to_scan) == 0:
            raise Exception("Could not recognise any peptides from input file")

    # list loaded peptides
    if verbose:
        print("Input epitopes:")
        for epitope in epitopes_to_scan:
            print(f"{epitope.name} {epitope.sequence}")
        print()

    return epitopes_to_scan

def MapPeptides(peptides, proteome, ref_genome, verbose=True):
    """
    Map peptides onto reference proteome and assign coding DNA sequences

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

        # if peptide was not prelocated from protein coordinates
        if peptide.genome_start is None:
            peptide.LocateSubsequence(proteome)

        # check if peptide was located successfully
        if peptide.genome_start is None:
            couldnt_map += 1
        # check ambiguous location
        elif peptide.genome_start == -1:
            mapped_ambiguously += 1
        # for mapped peptide
        else:
            # assign coding DNA sequence for mapped peptide
            peptide.AssignCodingSequence(ref_genome)
            mapped_peptides.append(peptide)

    if len(mapped_peptides) == 0:
        raise Exception('Could not map any peptides')

    # list mapped peptides and unmapped cases
    if verbose:
        print(f"Mapped {len(mapped_peptides)} peptides")
        print(f"Could not map {couldnt_map}")
        print(f"Mapped ambiguously {mapped_ambiguously}\n")

        for peptide in mapped_peptides:
            print(peptide)
        print()

    return mapped_peptides
