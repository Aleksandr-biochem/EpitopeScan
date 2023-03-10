#! /usr/bin/env python
import re
import os
import sys
import time
sys.path.append(os.path.realpath(os.path.dirname(__file__)))
from utils import *
import pandas as pd

if __name__ == "__main__":

	import argparse

	parser = argparse.ArgumentParser(description="EpitopeScan. Scan and analyse mutations within SARS-CoV-2 immunigenic peptides from Multiple Sequence Alignment")

	subparsers = parser.add_subparsers(title='Mode', dest='mode')
	
	# define scan mode options
	parser_scan = subparsers.add_parser('scan', help='Scan MSA file for epitope mutations')
	parser_scan.add_argument("-e", "--epitope", help="Input peptide epitope name and sequence comma-separated. E.x. S1,VGYWA", type=str)
	parser_scan.add_argument("-f", "--file",    help="Or path to input fasta file with multiple peptide epitopes", type=str)
	parser_scan.add_argument("-m", "--msa",     help="Path to MSA fasta file for analysis", type=str, required=True)
	parser_scan.add_argument("-o", "--out",     help="Output directory", type=str)
	parser_scan.add_argument("-t", "--tag",     help="Sample tag to filter", type=str)
	parser_scan.add_argument("-q", "--quality_filter", help="Threshold of max N bases proportion in genome. Recommended 0.05", type=float)
	parser_scan.add_argument("-n", "--no_ambiguity",   help="Discard any sample with ambiguous bases in epitope region", action='store_true')
	parser_scan.add_argument("-b", "--blosum",  help="Define BLOSUM version for mutation scoring. Default 90", type=int, default=90)
	parser_scan.add_argument("-s", "--sort",    help="Sort mutations in summary by count(0) or score(1). Default 0", type=int, choices=[0, 1], default=0)
	parser_scan.add_argument("-a", "--stat",    help="Stat individual mutations(0) or mutation combinations(1). Default 0", type=int, choices=[0, 1], default=0)
	
	# define stat mode options
	parser_stat = subparsers.add_parser('stat', help='Read existing output directory and print stats')	
	parser_stat.add_argument("-i", "--input",   help="Direcory with scan output to stat", type=str, required=True)
	parser_stat.add_argument("-b", "--blosum",  help="Define BLOSUM version for mutation scoring. Default 90", type=int, default=90)
	parser_stat.add_argument("-s", "--sort",    help="Sort mutations in summary by count(0) or score(1). Default 0", type=int, choices=[0, 1], default=0)
	parser_stat.add_argument("-a", "--stat",    help="Stat individual mutations(0) or mutation combinations(1). Default 0", type=int, choices=[0, 1], default=0)

	args = parser.parse_args()

	# load reference genome
	script_dir = os.path.realpath(os.path.dirname(__file__))
	with open(f'{script_dir}/reference_sequences/EPI_ISL_402124.fasta', 'r') as f:
		lines = f.readlines()
	reference_genome = lines[1].strip()
	del lines

	# load reference proteome
	proteome = ReadProteinsFromFile(f'{script_dir}/reference_sequences/protein_sequences_reference.fasta')

	# work in scan mode
	if args.mode == 'scan':

		# handle epitope inputs
		if not args.epitope is None:
			name, seq = args.epitope.split(',')
			if len(set(seq) - set('GALMFWKQESPVICYHRNDT')) > 0:
				raise Exception("Epitope sequence contains unrecognised characters")
			epitopes_to_scan = [Protein(name, seq)]
			print(f"Input epitope {epitopes_to_scan[0].name} {epitopes_to_scan[0].sequence}\n")
		elif not args.file is None:
			# read protein sequences from fasta file
			epitopes_to_scan = ReadProteinsFromFile(args.file)
			if len(epitopes_to_scan) > 0:
				print(f"Found {len(epitopes_to_scan)} input epitopes in {args.file}\n")
			else:
				raise Exception("Could not recognise any epitopes from input file")
		else:
			raise Exception("No epitopes provided")

	    # check if MSA file exists
		if not os.path.exists(args.msa):
			raise Exception(f"Provided MSA file does not exist, check path")

		# map epitopes onto proteome and assign coding DNA sequences
		epitopes_to_scan  = MapPeptides(epitopes_to_scan,
									 	proteome,
									  	reference_genome)

		# compile regex pattern to filter samples if any
		sample_tag = re.compile(args.tag) if args.tag else None

		# scan the MSA data
		print("Scanning MSA data for epitope mutations...")
		start_time = time.time()
		output_data = ScanMSA(epitopes_to_scan = epitopes_to_scan,
							  msa_file = args.msa,
							  sample_tag = sample_tag,
							  quality_filter = args.quality_filter,
							  ambiguity_intolerance = args.no_ambiguity)
		print(f"Finished scan. Run time: {round(time.time() - start_time, 2)} s.\n")

		# print key summary
		for epitope, data in zip(epitopes_to_scan, output_data):
			MutationTextSummary(epitope,
								data,
								stat_mutations=args.stat,
								sort_by=args.sort,
								blosum_version=args.blosum)

		# create output direcrories and save files
		output_dir = args.out if args.out else f"EpitopeScan_{time.strftime('%Y%m%d_%H%M%S')}"
		os.makedirs(output_dir, exist_ok=True)
		for epitope, data in zip(epitopes_to_scan, output_data):
			os.makedirs(f"{output_dir}/{epitope.name}", exist_ok=True) # create epitope subdir
			
			# save AA substitution matrix
			matrix = pd.DataFrame(epitope.AA_mutations_matrix,
					              index=list(epitope.sequence),
					              columns=list('GALMFWKQESPVICYHRNDT*'))
			matrix.to_csv(f"{output_dir}/{epitope.name}/{epitope.name}_AA_mutation_matrix.csv")

			# save NA substitution matrix
			matrix = pd.DataFrame(epitope.NA_mutations_matrix,
					              index=list(epitope.coding_sequence),
					              columns=list('ATCG'))
			matrix.to_csv(f"{output_dir}/{epitope.name}/{epitope.name}_NA_mutation_matrix.csv")
			
			# save output dataframes
			data.to_csv(f"{output_dir}/{epitope.name}/{epitope.name}_mutation_data.tsv", sep='\t')

		print(f"Saved outputs in {output_dir}")

	# work in stat mode
	if args.mode == 'stat':

		print(f"Collecting the data from {args.input}")
		
		StatEpitopeData(args.input,
						proteome,
						reference_genome,
						stat_mutations=args.stat,
						sort_by=args.sort,
						blosum_version=args.blosum)

