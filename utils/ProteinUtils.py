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






