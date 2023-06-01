#! /usr/bin/env python
import re
import os
import sys
import time
import numpy as np
import pandas as pd
from datetime import datetime

sys.path.append(os.path.realpath(os.path.dirname(__file__)))
from utils.Stat import *
from utils.SeqAnalysis import *
from utils.ProteinUtils import *

if __name__ == "__main__":

	import argparse

	parser = argparse.ArgumentParser(description="EpitopeScan. Scan and analyse SARS-CoV-2 genome Multiple Sequence Alignment for mutations within peptides")

	# create subparsers for two operation modes
	subparsers = parser.add_subparsers(title='Mode', dest='mode')
	
	# define scan mode key options
	parser_scan = subparsers.add_parser('scan', help='Scan genome MSA file for peptide mutations')
	parser_scan.add_argument("-e", "--epitope", help="Input peptide epitope name and sequence, comma-separated. E.x. S1,VGYWA", type=str)
	parser_scan.add_argument("-f", "--file",    help="Or path to fasta file with multiple input peptides", type=str)
	parser_scan.add_argument("--msa", 			help="Path to input MSA fasta file", type=str, required=True)
	parser_scan.add_argument("--metadata",    	help="PAth to metadata csv file to merge with mutaion data", type=str)
	
	# scan mode additional options
	parser_scan.add_argument("-o", "--out",     help="Output directory name", type=str)
	parser_scan.add_argument("-t", "--tag",     help="Sample tag to filter", type=str)
	parser_scan.add_argument("-q", "--quality_filter", help="Max threshold of N bases proportion in genome. Recommended 0.05", type=float)
	parser_scan.add_argument("-n", "--no_ambiguity",   help="Treat presence of any ambiguous bases in peptide region as insufficient coverage", action='store_true')
	parser_scan.add_argument("-b", "--blosum",  help="BLOSUM version for mutation scoring. Default 90", type=int, default=90)
	parser_scan.add_argument("-s", "--sort",    help="Sort mutations summary by count(0) or score(1). Default 0", type=int, choices=[0, 1], default=0)
	parser_scan.add_argument("-a", "--stat",    help="Stat individual mutations(0) or combinations(1). Default 0", type=int, choices=[0, 1], default=0)
	parser_scan.add_argument("--stat_with_metadata", help="Only stat samples with metadata", action='store_true')


	# define stat mode options
	parser_stat = subparsers.add_parser('stat', help='Read preexisting output and print mutation stats')	
	parser_stat.add_argument("-i", "--input",   help="Direcory with scan output", type=str, required=True)
	parser_stat.add_argument("-b", "--blosum",  help="BLOSUM version for mutation scoring. Default 90", type=int, default=90)
	parser_stat.add_argument("-s", "--sort",    help="Sort mutations summary by count(0) or score(1). Default 0", type=int, choices=[0, 1], default=0)
	parser_stat.add_argument("-a", "--stat",    help="Stat individual mutations(0) or combinations(1). Default 0", type=int, choices=[0, 1], default=0)
	parser_stat.add_argument("--stat_with_metadata", help="Only stat samples with metadata", action='store_true')
	parser_stat.add_argument("--start_date",    help="Stat samples after this date, dd/mm/yyyy", type=str)
	parser_stat.add_argument("--end_date",      help="Stat samples before this date, dd/mm/yyyy", type=str)

	args = parser.parse_args()

	## load reference genome 
	script_dir = os.path.realpath(os.path.dirname(__file__))
	with open(f'{script_dir}/reference_sequences/EPI_ISL_402124.fasta', 'r') as f:
		lines = f.readlines()
	reference_genome = lines[1].strip()
	del lines

	## load reference proteome
	proteome = ReadProteinsFromFile(f'{script_dir}/reference_sequences/protein_sequences_reference.fasta')

	## operation in scan mode
	if args.mode == 'scan':

		## check peptide inputs

		# single input peptide
		if not args.epitope is None:
			name, seq = args.epitope.split(',')
			if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
				raise Exception("Peptide sequence contains unrecognised characters")
			epitopes_to_scan = [Protein(name, seq)]
			print(f"Input epitope {epitopes_to_scan[0].name} {epitopes_to_scan[0].sequence}\n")

		# multiple peptides as fasta
		elif not args.file is None:
			epitopes_to_scan = ReadProteinsFromFile(args.file)
			if len(epitopes_to_scan) > 0:
				print(f"Found {len(epitopes_to_scan)} input peptides in {args.file}\n")
			else:
				raise Exception("Could not recognise any peptides from input file")

		# if no input epitopes provided
		else:
			raise Exception("No epitopes provided")

		# map epitopes onto proteome and assign coding DNA sequences
		epitopes_to_scan  = MapPeptides(epitopes_to_scan,
									 	proteome,
									  	reference_genome)

		# scan MSA data
		print("Scanning MSA data for mutations...")
		start_time = time.time()
		output_data = ScanMSA(epitopes_to_scan = epitopes_to_scan,
							  msa_file = args.msa,
							  sample_tag = sample_tag,
							  quality_filter = args.quality_filter,
							  ambiguity_intolerance = args.no_ambiguity)
		print(f"Finished scan. Run time: {round(time.time() - start_time, 2)} s.\n")

		## bind output DataFrames to MetaData
		BindMetadata(output_data, args.metadata, sample_tag)		

		## print key mutation summary
		if args.stat_with_metadata:
			print("Only calculating mutation stats for samples with metadata\n")

		for epitope, data in zip(epitopes_to_scan, output_data):
			df = data[data['has_metadata'] == 1] if args.stat_with_metadata else data
			summary = MutationTextSummary(epitope,
										  df,
										  stat_mutations=args.stat,
										  sort_by=args.sort,
										  blosum_version=args.blosum)
			print(f"{summary.shape[0]} mutations for {epitope.name}")
			if summary.shape[0] > 0:
				print(summary.to_string())
			print()

		## create output directories and save files
		print("Saving output...")
		# create main out dir
		if len(epitopes_to_scan) == 1:
			default_name = f"EpitopeScan_{epitopes_to_scan[0].name}_{time.strftime('%Y%m%d_%H%M%S')}"
		else:
			default_name = f"EpitopeScan_{time.strftime('%Y%m%d_%H%M%S')}"
		output_dir = args.out if args.out else default_name
		os.makedirs(output_dir, exist_ok=True)
		os.chdir(f"./{output_dir}")

		for epitope, data in zip(epitopes_to_scan, output_data):

			# for multiple epitopes create subdirs
			if len(epitopes_to_scan) > 1:
				os.makedirs(epitope.name, exist_ok=True) # create epitope subdir
				os.chdir(f"./{epitope.name}")

			# save AA substitution matrix
			matrix = pd.DataFrame(epitope.AA_mutations_matrix,
					              index=list(epitope.sequence),
					              columns=list('GALMFWKQESPVICYHRNDT*Δ'))
			matrix.to_csv(f"{epitope.name}_AA_mutation_matrix.csv")

			# save NA substitution matrix
			matrix = pd.DataFrame(epitope.NA_mutations_matrix,
					              index=list(epitope.coding_sequence),
					              columns=list('ATCGΔ'))
			matrix.to_csv(f"{epitope.name}_NA_mutation_matrix.csv")
			
			# save output mutation data
			data.to_csv(f"{epitope.name}_mutation_data.tsv", sep='\t')

			if len(epitopes_to_scan) > 1:
				os.chdir("../")

		print(f"Saved outputs to {output_dir}")

	## operation in stat mode
	elif args.mode == 'stat':

		print(f"Collecting mutation data from {args.input}...\n")

		start_date = datetime.strptime(args.start_date, '%d/%m/%Y') if not args.start_date is None else None
		end_date   = datetime.strptime(args.end_date, '%d/%m/%Y') if not args.end_date is None else None

		StatEpitopeData(args.input,
						proteome,
						reference_genome,
						stat_mutations=args.stat,
						sort_by=args.sort,
						time_start=start_date,
						time_end=end_date,
						blosum_version=args.blosum,
						metadata_filter=args.stat_with_metadata)

