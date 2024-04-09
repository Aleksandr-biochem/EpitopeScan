import sys
import os
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from EpitopeScan.utils.ProteinUtils import *

def VerifyInput(files):
	"""
	Verify input files for EpitopeScanGUI

	files - streamlit file_uploader
	
	Returns:
	right_input - bool, True if inputs are right
	""" 

	right_input = False

	if len(files) > 3: 
		st.sidebar.write(f"Expected 3 files, but {len(files)} provided. Delete extra")
	elif len(files) < 3: 
		st.sidebar.write(f"Expected 3 files, but {len(files)} provided. Upload missing")
	else:
		# get epitope neames from files
		file_names = [files[i].name for i in range(len(files))]
		names = set([name.split('_')[0] for name in file_names])
		if len(names) > 1:
			st.sidebar.write(f"Input files seem to belong to different epitopes. Revise input")
		else:
			name = names.pop()
			if not f"{name}_AA_mutation_matrix.csv" in file_names:
				st.sidebar.write(f"Missing {name}_AA_mutation_matrix.csv")
			elif not f"{name}_NA_mutation_matrix.csv" in file_names:
				st.sidebar.write(f"Missing {name}_NA_mutation_matrix.csv")
			elif not f"{name}_mutation_data.tsv" in file_names:
				st.sidebar.write(f"Missing {name}_mutation_data.tsv")
			else:
				right_input = True

	return right_input

@st.cache_data
def LoadData(files):
	"""
	Load EpitopeScan outputs into DataFrames

	files - streamli file_uploader

	Returns:
	peptide_name, mutation_data, AA_mutation_matrix, NA_mutation_matrix
	"""
	
	# get peptide name
	peptide_name = files[0].name.split('_')[0]

	for i in range(len(files)):

		if files[i].name == f"{peptide_name}_AA_mutation_matrix.csv":
			AA_mutation_matrix = pd.read_csv(files[i], index_col=0)
			peptide = Protein(peptide_name, ''.join(AA_mutation_matrix.index)) 

		elif files[i].name == f"{peptide_name}_NA_mutation_matrix.csv":
			NA_mutation_matrix = pd.read_csv(files[i], index_col=0)

		elif files[i].name == f"{peptide_name}_mutation_data.tsv":
			mutation_data = pd.read_csv(files[i], sep='\t', index_col=0)
			mutation_data['sample_date']= pd.to_datetime(mutation_data['sample_date'],
														 format='%Y-%m-%d')

	return peptide, mutation_data, AA_mutation_matrix, NA_mutation_matrix

@st.cache_data
def KeyPieChart(key_stat):
	"""
	Assemble pie chart of key sample statistics

	key_stat - dict, key samples statistics

	Returns:
	fig - go.Figure, assembled figure
	"""
	# create figure
	fig = go.Figure()

	fig.add_trace(
		go.Pie(labels=list(key_stat.keys()),
			   values=list(key_stat.values()))
				 )

	fig.update_traces(textfont_size=20)
	fig.update_layout(title_text = "Samples breakdown",
					  title_x = 0.5,
					  title_font_size=20,
					  legend_font_size=15,
					  height=500, width=500)

	return fig

@st.cache_data
def KeyWeeklyStat(df, report_proportion=False):
	"""
	Calculate key statistics per week

	df - pd.DataFrame, mutation data to aggregate
	report_proportion - bool, report proportion instead of counts
	
	Returns:
	stat - pd.DataFrame, key statistics per week
	"""

	# create dummy week column
	max_week = df['epi_week'].max()
	min_week = df['epi_week'].min()
	week_column = pd.DataFrame()
	week_column['epi_week'] = [w for w in range(min_week, max_week + 1)]

	# get statistics for each type of reports

	# non-functional
	grouped_nf = df[df['AA_mutations'] == "NF"] \
	             .groupby('epi_week', as_index=False) \
	             .count()[["epi_week", 'sequence_name']] \
	             .rename(columns={'sequence_name':'Non-functional'})

	# no coverage
	grouped_nc = df[df['AA_mutations'] == "-"] \
	             .groupby('epi_week', as_index=False) \
	             .count()[["epi_week", 'sequence_name']] \
	             .rename(columns={'sequence_name':'No coverage'})

	# no mutations
	grouped_nm = df[df['AA_mutations'].isna()] \
	             .groupby('epi_week', as_index=False) \
	             .count()[["epi_week", 'sequence_name']] \
	             .rename(columns={'sequence_name':'No mutations'})

	# with mutations
	grouped_mu = df[(df['AA_mutations'] != "-")  & \
	             (df['AA_mutations'] != "NF") & \
	             (~df['AA_mutations'].isna())] \
	             .groupby('epi_week', as_index=False) \
	             .count()[["epi_week", 'sequence_name']] \
	             .rename(columns={'sequence_name':'With mutations'})

	# total samples
	grouped_to = df.groupby('epi_week', as_index=False) \
	             .count()[["epi_week", 'sequence_name']] \
	             .rename(columns={'sequence_name':'Total samples'})

	# merge all DataFrames
	stat = week_column.merge(grouped_nf, on='epi_week', how='left') \
	       .merge(grouped_nc, on='epi_week', how='left') \
	       .merge(grouped_nm, on='epi_week', how='left') \
	       .merge(grouped_mu, on='epi_week', how='left') \
	       .merge(grouped_to, on='epi_week', how='left') \

	stat = stat.fillna(0).astype(int)

	if report_proportion:
	    for col in stat.columns[1:-1]:
	        stat[col] = stat[col] / stat['Total samples']
	    stat.drop(columns=['Total samples'], inplace=True)

	return stat

@st.cache_data
def PlotKeyStat(mutation_data, plot_proportion=False):
	"""
	Assemble bar chart of sample statistics in time

	mutation_data - pd DataFrame, mutation data for plotting
	plot_proportion - bool, plot as proportion of all samples

	Returns:
	fig - go.Figure, assembled figure
	"""
	mutation_data = mutation_data[mutation_data['has_metadata'] == 1]
	mutation_data['epi_week'] = mutation_data['epi_week'].astype('int32')

	# get statistics
	stat = KeyWeeklyStat(mutation_data, plot_proportion)

	# create figure
	fig = go.Figure()

	# bar plots for each statistic
	for value in stat.columns[1:]:
	    fig.add_trace(go.Bar(x=stat['epi_week'],
							 y=stat[value],
							 name=value)
					 )

	fig.update_layout(title_text = "Samples statistics per week",
					  title_x = 0.5,
					  xaxis_title = "Pandemic week",
					  height=500, width=800)

	if plot_proportion:
	    fig.update_layout(yaxis_title = "Proportion in total samples per week")
	else:
	    fig.update_layout(yaxis_title = "Count")
	    
	return fig

@st.cache_data
def StatMutations(df,
                  mutations,
                  report_proportion=False,
                  time_start=None,
                  time_end=None):
    """
    Calculate statistics per mutation

    df - pd DataFrame, mutation data to stat
    mutations - list(str), list of mutations to stat
    report_proportion - bool, report proportion instead of counts
    time_start - datetime, start date to filter
    time_end - datetime, end time to filter

    Returns:
    stat_df - pd DataFrame, calculated statistic per mutation
    """
    # convert week to int
    df['epi_week'] = df['epi_week'].astype('int32')
    
    # filter by date if necessary
    if not time_start is None:
        df = df[df['sample_date'] >= time_start]
    if not time_end is None:
        df = df[df['sample_date'] <= time_end]
        
    # create dummy week column
    max_week = df['epi_week'].max()
    min_week = df['epi_week'].min()
    week_column = pd.DataFrame()
    week_column['epi_week'] = [w for w in range(min_week, max_week + 1)]
        
    # drop samples without coverage
    df = df[df['AA_mutations'] != "-"]
    
    # count samples per week
    sample_counts = df.groupby('epi_week', as_index=False) \
                    .count()[["epi_week", 'sequence_name']] \
                    .rename(columns={'sequence_name': 'n_samples'})
    
    
    # get statistics for each mutation (combination)
    counts_per_mutation = []
    for mut in mutations:
        
        muts = mut.split(',')
        
        # subset samples with target mutation
        mask = ~df['sequence_name'].isna()
        for m in muts:
            mask = mask & (df[m] == 1)

        # impose additional filtering for exact combination
        if len(muts) > 1:
            unwanted_muts = set(df.columns[3:-4]) - set(muts)
            for m in unwanted_muts:
                   mask = mask & (df[m] == 0)

        # count samples with target mutation/combination
        mutation_counts = df[mask].groupby('epi_week', as_index=False) \
                          .count()[['epi_week', 'sequence_name']] \
                          .rename(columns={'sequence_name': mut})
        
        counts_per_mutation.append(mutation_counts)
    
    # merge all DataFrames
    stat_df = week_column.merge(sample_counts, on='epi_week', how='left')
    for counts in counts_per_mutation:
        stat_df = stat_df.merge(counts, on='epi_week', how='left')
    
    stat_df = stat_df.fillna(0).astype(int)
    
    if report_proportion:
        for col in stat_df.columns[2:]:
            stat_df[col] = stat_df[col] / stat_df['n_samples']
        stat_df.drop(columns=['n_samples'], inplace=True)
    
    return stat_df

@st.cache_data
def PlotMutationStat(df,
                     mutations,
                     plot_proportion=False,
                     time_start=None, time_end=None):
    """
    Assemble bar chart of mutation statistics per week

    df - pd DataFrame, mutation data to stat
    mutations - list(str), list of mutations to stat
    report_proportion - bool, report proportion instead of counts
    time_start - datetime, start date to filter
    time_end - datetime, end time to filter

    Returns:
    fig - go.Figure, assembled figure
    """
    
    # get statistics
    stat = StatMutations(df, mutations, plot_proportion)

    # create figure
    fig = go.Figure()

    # bar plot for mutation statistic
    i = 1 if plot_proportion else 2
    for col in stat.columns[i:]:
        fig.add_trace(go.Bar(x=stat['epi_week'],
                             y=stat[col],
                             name=col))
    
    fig.update_layout(title_text = f"Mutations statistics per week",
                      title_x = 0.5,
                      xaxis_title = "Pandemic week")

    if plot_proportion:
        fig.update_layout(yaxis_title = "Proportion in samples with coverage")
    else:
        fig.add_trace(go.Bar(x=stat['epi_week'],
                         y=stat['n_samples'],
                         name='Total samples'))
        fig.update_layout(yaxis_title = "Count")

    return fig

@st.cache_data
def PlotSeqConservation(NA_matrix, AA_matrix):
    """
    Assemble bar chart of overall AA and NA sequence conservation

    NA_matrix - pd DataFrame, NA substitution count matrix
    AA_matrix - pd DataFrame, AA substitution count matrix

    Returns:
    fig - go.Figure, assembled figure
    """
    
    # calculate conservation
    NA_matrix = NA_matrix.apply(lambda row: row/row.sum(), axis=1)
    AA_matrix['Conservation'] = AA_matrix.apply(lambda row: row[row.name]/row.sum(), axis=1)
    
    # plot conservation
    fig = make_subplots(rows=2, cols=1,
                        y_title = "Conservation",
                        x_title = "Sequence",
                        vertical_spacing=0.08)

    # AA conservation
    fig.append_trace(go.Bar(x = [i for i in range(len(AA_matrix.index))],
                            y = AA_matrix['Conservation'],
                            name="Protein sequence"),
                     row=1, col=1)

    fig.update_xaxes(tickmode = 'array',
                     tickvals = [i for i in range(len(AA_matrix.index))],
                     ticktext = AA_matrix.index.values,
                     row=1, col=1)

    # NA conservation
    for b in 'ATCGΔ':
        fig.append_trace(go.Bar(x = [i for i in range(len(NA_matrix.index))] ,
                                y = NA_matrix[b],
                                name = b),
                         row=2, col=1)

    fig.update_xaxes(tickmode = 'array',
                     tickvals = [i for i in range(len(NA_matrix.index))],
                     ticktext = NA_matrix.index.values,
                     row=2, col=1)

    fig.update_layout(title_text="Sequence conservation",
                      title_x=0.5,
                      showlegend=False,
                      barmode='stack',
                      height=500, width=800)

    return fig

@st.cache_data
def StatLineages(df,
                 mutations,
                 report_proportion,
                 time_start=None, time_end=None):
    """
    Aggreagate samples data by lineage and mutations

    df - pd DataFrame, mutation data to stat
    mutations - list(str), list of mutations to stat
    report_proportion - bool, report proportion instead of counts
    time_start - datetime, start date to filter
    time_end - datetime, end time to filter

    Returns:
    lineage_counts - pd.DataFrame, sample data aggregated by lineage and mutation
    """
    
    # filter by date if necessary
    if not time_start is None:
        df = df[df['sample_date'] >= time_start]
    if not time_end is None:
        df = df[df['sample_date'] <= time_end]
    
    # count lineages present in this subset
    lineage_counts = df.groupby('usher_lineage', as_index=False) \
                      .count()[["usher_lineage", 'sequence_name']] \
                      .rename(columns={'sequence_name': 'n_samples'})
    

    # get statistics for all mutations
    counts_per_mutation = []
    for mut in mutations:
        
        muts = mut.split(',')
        
        # subset samples with target mutation
        mask = ~df['sequence_name'].isna()
        for m in muts:
            mask = mask & (df[m] == 1)

        # impose additional filtering for exact combination
        if len(muts) > 1:
            unwanted_muts = set(df.columns[3:-4]) - set(muts)
            for m in unwanted_muts:
                   mask = mask & (df[m] == 0)

        # count samples with target mutations
        mutation_counts = df[mask].groupby('usher_lineage', as_index=False) \
                          .count()[["usher_lineage", 'sequence_name']] \
                          .rename(columns={'sequence_name': mut})
        
        if mutation_counts[mut].sum() > 0:
        	counts_per_mutation.append(mutation_counts)
    
    # merge all DataFrames
    for counts in counts_per_mutation:
        lineage_counts = lineage_counts.merge(counts, on='usher_lineage', how='left')
    
    lineage_counts = lineage_counts.fillna(0)
    for col in lineage_counts.columns[1:]:
        lineage_counts[col] = lineage_counts[col].astype(int)
    
    if report_proportion:
        for col in lineage_counts.columns[2:]:
            lineage_counts[col] = lineage_counts[col] / lineage_counts['n_samples']
        lineage_counts.drop(columns=['n_samples'], inplace=True)
    
    return lineage_counts