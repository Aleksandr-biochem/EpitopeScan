#! /usr/bin/env python
import os
import json
import subprocess
import shutil
import pandas as pd
import numpy as np

"""
run EpitopeScan tests on sample peptide sequences
"""
def LoadOutput(out_dir):

	name = out_dir.split('_')[-1]

	mut_data  = pd.read_csv(f"{out_dir}/{name}_mutation_data.tsv", sep='\t', index_col=0)
	NA_matrix = pd.read_csv(f"{out_dir}/{name}_NA_mutation_matrix.csv", index_col=0)
	AA_matrix = pd.read_csv(f"{out_dir}/{name}_AA_mutation_matrix.csv", index_col=0)

	return mut_data, NA_matrix, AA_matrix

def CompareMutationData(mut_df, ref_data, ambiguity_intolerance=False):
	
	mismatched_samples = []

	for _, row in mut_df.iterrows():
		
		name = row['sequence_name']

		if row['NA_mutations'] == 'NF':
			if not ref_data[name]['non_functional']:
				mismatched_samples.append(name)
				print(f"{name} reported as NF, when it is not")

		elif row['NA_mutations'] == '-':
			if not ref_data[name]['no_coverage']:
				if ambiguity_intolerance and ref_data[name]['has_ambiguity']:
					continue
				else: 
					mismatched_samples.append(name)
					print(f"{name} reported as no coverage, when it should not")
		else:

			if ref_data[name]['non_functional']:
				mismatched_samples.append(name)
				print(f"{name} not reported as NF")

			elif ref_data[name]['no_coverage']:
				mismatched_samples.append(name)
				print(ref_data[name]['no_coverage'])
				print(f"{name} not reported as no coverage")

			elif row['NA_mutations'] is np.nan:
				if not len(ref_data[name]['NA_mutations']) == 0:
					mismatched_samples.append(name)
					print(f"{name} has mutations when none expected")

			elif row['AA_mutations'] is np.nan:
				if not len(ref_data[name]['AA_mutations']) == 0:
					mismatched_samples.append(name)
					print(f"{name} has mutations when none expected")

			else:
				NA_mutations = set(row['NA_mutations'].split(','))
				AA_mutations = set(row['AA_mutations'].split(','))
				if NA_mutations != set(ref_data[name]['NA_mutations']):
					print(f"{name} NA mutation mismatch:", NA_mutations, ' ', ref_data[name]['NA_mutations'])
					mismatched_samples.append(name)
				if AA_mutations != set(ref_data[name]['AA_mutations']):
					print(f"{name} AA mutation mismatch:", AA_mutations, ' ', ref_data[name]['AA_mutations'])
					mismatched_samples.append(name)

	return mismatched_samples

if __name__ == "__main__":

	#### peptides to run tests on #####
	test_peptides = [
		('T1', 'IIRGWIFGTTLDSKTQSLLIV'),
		('T2', 'QSFLNRVCGVSAARL'),
		('T3', 'LTSHTVMPLSAPTLVPQEHY')
	]
	###################################

	# run tests for each peptide
	for name, peptide in test_peptides:

		print(f"Running tests on peptide {name} {peptide}\n")

		# load reference data
		script_dir = os.path.realpath(os.path.dirname(__file__))
		with open(f'{script_dir}/{name}_ref/{name}_data.json') as json_file:
			ref_data = json.load(json_file)
		ref_NA_matrix = pd.read_csv(f"{script_dir}/{name}_ref/{name}_NA_mutation_matrix.csv", index_col=0)
		ref_AA_matrix = pd.read_csv(f"{script_dir}/{name}_ref/{name}_AA_mutation_matrix.csv", index_col=0)

		
		### TEST 1 ###
		# default run options

		print(f"Running test 1...")

		# run analysis at default parameters
		subprocess.run([f"{script_dir}/../EpitopeScan.py", "scan",
						"-e", f"{name},{peptide}",
						"--msa", f"{script_dir}/{name}_ref/{name}_test_sequences.fasta",
						"-o", f"ES_{name}"],
						stdout=subprocess.DEVNULL)


		# load output
		mut_data, NA_matrix, AA_matrix = LoadOutput(f"ES_{name}")

		# compare count matrices
		NA_count_comparison = NA_matrix.compare(ref_NA_matrix)
		AA_count_comparison = AA_matrix.compare(ref_AA_matrix)

		# compare mutation matrices
		mistmatches = CompareMutationData(mut_data, ref_data)
		if (len(mistmatches) == 0) and NA_count_comparison.empty and AA_count_comparison.empty:
			print("Test passed\n")

		############

		### TEST 2 ###
		# tag filtering

		print(f"Running test 2...")

		# run analysis at default parameters with filtering
		subprocess.run([f"{script_dir}/../EpitopeScan.py", "scan",
						"-e", f"{name},{peptide}",
						"--msa", f"{script_dir}/{name}_ref/{name}_test_sequences.fasta",
						"-t", "_tag",
						"-o", f"ES_{name}"],
						stdout=subprocess.DEVNULL)


		# load output
		mut_data, NA_matrix, AA_matrix = LoadOutput(f"ES_{name}")

		# compare mutation matrices
		mistmatches = CompareMutationData(mut_data, ref_data)

		# check that no samples to be filtered out passed
		samples_which_should_be_filtered = sum(mut_data['sequence_name'].str.match('out'))
		samples_which_should_be_kept = sum(~mut_data['sequence_name'].str.match('out'))
		samples_which_should_be_kept_ref = sum([1 if 'tag' in name else 0 for name in ref_data])
		if (len(mistmatches) == 0) and \
		   (samples_which_should_be_filtered == 0) and \
		   (samples_which_should_be_kept == samples_which_should_be_kept_ref):
			print("Test passed\n")
		############

		### TEST 3 ###
		# quality filtering

		print(f"Running test 3...")

		# run analysis with quality filter
		subprocess.run([f"{script_dir}/../EpitopeScan.py", "scan",
						"-e", f"{name},{peptide}",
						"--msa", f"{script_dir}/{name}_ref/{name}_test_sequences.fasta",
						"-q", "0.05",
						"-o", f"ES_{name}"],
						stdout=subprocess.DEVNULL)


		# load output
		mut_data, NA_matrix, AA_matrix = LoadOutput(f"ES_{name}")

		# check that all samples are correctly filtered
		should_be_kept1 = set()
		for seq_name in ref_data:
			if not ref_data[seq_name]['high_N_content']:
				should_be_kept1.add(seq_name)
		processed_samples1 = set(mut_data['sequence_name'])

		# this should return all samples
		subprocess.run([f"{script_dir}/../EpitopeScan.py", "scan",
						"-e", f"{name},{peptide}",
						"--msa", f"{script_dir}/{name}_ref/{name}_test_sequences.fasta",
						"-q", "0.30",
						"-o", f"ES_{name}"],
						stdout=subprocess.DEVNULL)


		# load output
		mut_data, NA_matrix, AA_matrix = LoadOutput(f"ES_{name}")

		# check that all samples are correctly filtered
		should_be_kept2 = set()
		for seq_name in ref_data:
			should_be_kept2.add(seq_name)
		processed_samples2 = set(mut_data['sequence_name'])

		if (processed_samples1 == should_be_kept1) and \
		   (processed_samples2 == should_be_kept2):
			print('Test passed\n')
		############

		### TEST 4 ###
		# ambiguity intolerance 
		
		print(f"Running test 4...")

		# run analysis with ambiguity intolerance
		subprocess.run([f"{script_dir}/../EpitopeScan.py", "scan",
						"-e", f"{name},{peptide}",
						"--msa", f"{script_dir}/{name}_ref/{name}_test_sequences.fasta",
						"-n",
						"-o", f"ES_{name}"],
						stdout=subprocess.DEVNULL)


		# load output
		mut_data, NA_matrix, AA_matrix = LoadOutput(f"ES_{name}")

		# compare mutation matrices
		mistmatches = CompareMutationData(mut_data, ref_data, ambiguity_intolerance=True)

		if len(mistmatches) == 0:
			print("Test passed\n")
		############

		# delete output
		shutil.rmtree(f'ES_{name}')












