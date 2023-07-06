#! /usr/bin/env python
import os
import sys
import pandas as pd
import streamlit as st
from datetime import datetime, timedelta 

sys.path.append(os.path.realpath(os.path.dirname(__file__)))
from utils.GUIutils import *
from utils.Stat import *
from utils.ProteinUtils import *

def main():
 
	# Headline
	st.markdown("""
	### EpitopeScan mutation data analysis

- This app facilitates the analysis of mutation data generated with `EpitopeScan`

- Built with `Python`, `Streamlit` and `Plotly`
---""")

	# Sidebar for input
	with st.sidebar.header('Input data'):
	
		st.sidebar.markdown("""
		Provide EpitopeScan outputs for analysis:
		""")

		# input files for analysis
		uploaded_files = st.sidebar.file_uploader("Upload 3 input files for epitope",
												  accept_multiple_files=True)

		# verify inputs
		verified = False
		if len(uploaded_files) > 0:
			verified = VerifyInput(uploaded_files)

	# proceed if input files are verified
	if verified:

		# load data from files
		with st.spinner(f"Loading input tables..."):
			peptide, mutation_data, AA_mutation_matrix, NA_mutation_matrix = LoadData(uploaded_files)

		# check if any samples have metadata
		no_metadata = True if sum(mutation_data['has_metadata'] == 1) == 0 else False

		## Print key sample statistics ##
		KeyStatisticsSection(peptide, mutation_data)

		## Plot key statistics in time ##
		KeyStatisticsPlotSection(mutation_data, no_metadata)

		## List and plot found mutations ##
		MutationStatSection(peptide, mutation_data, no_metadata)

		## Sequence conservation ##
		SequenceConservationSection(NA_mutation_matrix, AA_mutation_matrix)

		## Lineage statistics ##
		LineageStatSection(peptide, mutation_data, no_metadata)

def KeyStatisticsSection(peptide, mutation_data):
	"""
	App section with general sample statistics
	"""

	st.write('<p style="font-size:22px;weight:bold">1. General sample statistics</p>',
			 unsafe_allow_html=True)

	# load reference genome
	script_dir = os.path.realpath(os.path.dirname(__file__))
	with open(f'{script_dir}/reference_sequences/EPI_ISL_402124.fasta', 'r') as f:
		lines = f.readlines()
	reference_genome = lines[1].strip()
	del lines

	## load reference proteome
	proteome = ReadProteinsFromFile(f'{script_dir}/reference_sequences/protein_sequences_reference.fasta')

	# map peptide onto proteome and assign coding DNA sequence
	peptide = MapPeptides([peptide],
	                      proteome,
	                      reference_genome,
	                      verbose=False)[0]

	# print peptide info
	st.write(f"Peptide {peptide.name}, sequence: **{peptide.sequence}** ({len(peptide)} residues)")
	part1 = f"Originates from {peptide.parent_protein[0]} protein (residues {peptide.protein_start} to {peptide.protein_start + len(peptide) - 1})."
	part2 = f" Genome coordinates: {peptide.genome_start}-{peptide.genome_end}"
	st.write(part1 + part2)
	if len(peptide.parent_protein) > 1:
		st.write(f"Other parent proteins: {', '.join(peptide.parent_protein[1:])}")

	with st.container():
		# create two column layout
		col1, col2 = st.columns(2)

		# first column is to print statistics
		with col1:
			# whether to stat samples without metadata
			stat_only_metadata = st.checkbox('Only stat samples with metadata')
		
			key_stat = GetKeyStat(mutation_data,
				                  stat_only_metadata)

			# print statistics
			st.write(f"**Data for peptide {peptide.name} contains:**")
			st.write(f"{key_stat['Total']} samples")
			st.write(f"of them {key_stat['Has_metadata']} have metadata ({round(key_stat['Has_metadata']*100/key_stat['Total'], 1)}% of all)")
			
			# in case we can provide clear date range
			some_samples_lack_metadata = True if 0 in mutation_data['has_metadata'] else False
			if not some_samples_lack_metadata:
				start_date = mutation_data['sample_date'].min().strftime('%d/%m/%Y')
				end_date = mutation_data['sample_date'].max().strftime('%d/%m/%Y')
				st.write(f"dated between {start_date} and {end_date}")

			if stat_only_metadata:
				st.write("**Of samples with metadata:**")
			else:
				st.write("**Of all samples:**")

			total = key_stat['Has_metadata'] if stat_only_metadata else key_stat['Total']
			st.markdown(f"- {key_stat['Have_mutation']} ({round(key_stat['Have_mutation'] * 100/ total, 2)}%) have mutations")
			st.markdown(f"- {key_stat['No_mutation']} ({round(key_stat['No_mutation'] * 100/ total, 2)}%) have no mutations")
			st.markdown(f"- {key_stat['No_coverage']} ({round(key_stat['No_coverage'] * 100/ total, 2)}%) reported with insufficient coverage")
			st.markdown(f"- {key_stat['Non_functional']} ({round(key_stat['Non_functional'] * 100/ total, 2)}%) reported as non-functional")

		# second column to plot pie chart
		with col2:
			key_stat.pop('Has_metadata', None)
			key_stat.pop('Total', None)
			pie_chart = KeyPieChart(key_stat)
			st.plotly_chart(pie_chart)

def KeyStatisticsPlotSection(mutation_data, no_metadata):
	"""
	Plot sample categories in time
	as counts or proportions in weekly sample counts
	"""
	st.write('<p style="font-size:22px;weight:bold">2. General sample statistics in time</p>',
			 unsafe_allow_html=True)

	if no_metadata:
		st.write('<p style="font-size:18px">No samples with metadata to plot</p>',
			 	 unsafe_allow_html=True)
	else:
		# plot counts or proportion
		col1, col2 = st.columns((1, 3))
		with col1:
			plot_as2 = st.selectbox('Plot as:', ('counts', 'proportion'), key=1)
			plot_proportion2 = True if plot_as2 == 'proportion' else False

		with st.spinner(f"Generating weekly stat plot..."):
			key_stat_plot = PlotKeyStat(mutation_data,
										plot_proportion=plot_proportion2)
			st.plotly_chart(key_stat_plot)

def MutationStatSection(peptide, mutation_data, no_metadata):
	"""
	Section with statistics for each mutation
	"""
	st.write('<p style="font-size:22px;weight:bold">3. Stat mutations</p>',
			 unsafe_allow_html=True)

	# whether to stat samples without metadata
	stat_only_metadata = st.checkbox('Only stat samples with metadata',
									  value=True if not no_metadata else False,
									  key='stat_only_metadata')
	
	# additional parameters for stat calculation
	with st.container():
		col1, col2, col3, col4, col5 = st.columns((1.5, 1.5, 1.2, 1.2, 1))

		# check if date span specification is valid
		if not no_metadata:
			def_start_date = mutation_data['sample_date'].min()
			def_end_date = mutation_data['sample_date'].max()
		else:
			def_start_date, def_end_date = None, None

		with col1:
			stat_as = st.selectbox('Stat mutations as', ('individual', 'combinations'), key=2)
			stat_as_key = 0 if stat_as == 'individual' else 1
		with col2:
			blosum_v = st.text_input('BLOSUM score version', '90')
			blosum_v = int(blosum_v)
		with col3:
			start_date = st.date_input("Stat starting from:", def_start_date, key='d1',
										disabled=not stat_only_metadata)
			start_date = datetime(start_date.year, start_date.month, start_date.day)
		with col4:
			end_date = st.date_input("and up to:", def_end_date, key='d2',
									 disabled=not stat_only_metadata)
			end_date = datetime(end_date.year, end_date.month, end_date.day)
		with col5:
			plot_as3 = st.selectbox('Plot as:', ('counts', 'proportion'), key=3,
									disabled=not stat_only_metadata)
			plot_proportion3 = True if plot_as3 == 'proportion' else False

	# generate and display stats
	with st.container():

		# create two column layout
		col1, col2 = st.columns((1.7, 1))

		# print table mutation summary
		with st.spinner(f"Generating mutation summary..."):

			# filter mutation data by metadata and dates
			if stat_only_metadata:
				filtered_mutation_data = mutation_data[mutation_data['has_metadata'] == 1]
				filtered_mutation_data = filtered_mutation_data[filtered_mutation_data['sample_date'] >= start_date]
				filtered_mutation_data = filtered_mutation_data[filtered_mutation_data['sample_date'] <= end_date]
			else:
				filtered_mutation_data = mutation_data

			summary = MutationTextSummary(peptide,
										  filtered_mutation_data,
					                      stat_mutations=stat_as_key,
					                      sort_by=0,
					                      blosum_version=blosum_v)
			
			with col1:
				st.write(f"**Mutations summarised from {filtered_mutation_data.shape[0]} samples:**")
				st.write(summary)

		# mutation selection for plotting
		if stat_only_metadata:
				
				# select mutations to plot
				with col2:
					found_mutations = summary.iloc[:, 1].values
					top_found_mutations = found_mutations[:5] if len(found_mutations)>5 else found_mutations
					selected_mutations = st.multiselect('Mutations to plot:',
														found_mutations,
														top_found_mutations)
	# plot selected mutations in time
	if stat_only_metadata:
		with st.spinner(f"Generating timeline plot..."):
			# generate plot
			muts_time_plot = PlotMutationStat(mutation_data[mutation_data['has_metadata'] == 1],
						                      selected_mutations,
											  plot_proportion=plot_proportion3,
											  time_start=start_date,
											  time_end=end_date)
			
			st.plotly_chart(muts_time_plot, use_container_width=True)

@st.cache_data
def SequenceConservationSection(NA_mutation_matrix, AA_mutation_matrix):
	"""
	Plot sequence conservation
	"""
	with st.spinner(f"Generating sequence conservation plot..."):
		sequence_conservation_plot = PlotSeqConservation(NA_mutation_matrix,
														 AA_mutation_matrix)
		st.plotly_chart(sequence_conservation_plot, use_container_width=True)

def LineageStatSection(peptide, mutation_data, no_metadata):
	"""
	Generate statistics for samples lineages
	"""
	st.write('<p style="font-size:22px;weight:bold">4. Stat mutations by lineage</p>',
			 unsafe_allow_html=True)

	if no_metadata:
		st.write('<p style="font-size:18px">No samples with metadata to stat</p>',
			 	 unsafe_allow_html=True)
	else:
		# specify date range for analysis
		# three last weeks by default
		with st.container():
				col1, col2, col3, pad = st.columns((1.2, 1.2, 1, 1))

				def_end_date = mutation_data['sample_date'].max()
				def_start_date = def_end_date - timedelta(weeks=3) 

				with col1:
					start_date = st.date_input("Stat starting from:", def_start_date, key='d4')
					start_date = datetime(start_date.year, start_date.month, start_date.day)
				with col2:
					end_date = st.date_input("and up to:", def_end_date, key='d5')
					end_date = datetime(end_date.year, end_date.month, end_date.day)
				with col3:
					report_as = st.selectbox('Report samples as:', ('counts', 'proportion'), key=4)
					report_proportion = True if report_as == 'proportion' else False

		# select mutations to stat
		with st.spinner(f"Generating lineage statistics..."):

			summary2 = MutationTextSummary(peptide,
										   mutation_data[mutation_data['has_metadata'] == 1],
					                       stat_mutations=0,
					                       sort_by=0,
					                       blosum_version=90)
			found_mutations2 = summary2.iloc[:, 1].values

			# get lineage summary
			lineage_summary = StatLineages(mutation_data[mutation_data['has_metadata'] == 1],
							               mutations=found_mutations2,
							               report_proportion=report_proportion,
							               time_start=start_date,
							               time_end=end_date)
			st.write(lineage_summary)

if __name__ == "__main__":
	main()