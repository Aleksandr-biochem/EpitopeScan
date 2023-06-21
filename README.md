# EpitopeScan
EpitopeScan is a Python3 toolset facilitating the tracking of mutations in SARS-CoV-2 immunogenic epitopes.

**Contents:**
1. [Repository contents](#sec1) </br>
2. [Download and set-up](#sec2) </br>
3. [EpitopeScan usage](#sec3) </br>
    1. [Input data](#sec31) </br>
    1. [Running command line tool](#sec32) </br>
        1. [Scan mode](#sec321) </br>
        1. [Stat mode](#sec322) </br>
        1. [Mutation data output](#sec323) </br>
    1. [Running graphical user-interface tool](#sec33) </br>
    1. [Running test scripts](#sec34) </br>
4. [Algorithm details](#sec4) </br>

<a name="sec1"></a>
## 1. Repository contents

- `EpitopeScan.py` is a command-line tool for mutation report from SARS-CoV-2 multiple genome alignment and for calculation of general mutation statistics.

- `EpitopeScanGUI.py` is a Graphical User Interface (GUI) application. It uses system's default web-brwoser to provide interactive interface for the analysis of mutation data generated by `EpitopeScan.py`.

- `utils` folder contains supporting Python3 code with functions and classes for command-line and GUI tools.

- `reference_sequences` contains reference SARS-CoV-2 data listed from [GISAID](https://gisaid.org/wiv04/). This consists of reference genome `EPI_ISL_402124.fasta`, table with Open Reading Frames information `ORFs_reference.txt` (ORF name, genome start and end coordinates, name of translated protein, length of translated protein), reference protein sequences `protein_sequences_reference.fasta`.

- `test` contains subfolder with example data for analysis `example_data`, subfolders with test data `T#_...` and the script `run_EpitopeScan_test.py` for the command-line tool testing.

- `requirements.txt` is a list of packages to be installed when setting up Python3 environment.

<a name="sec2"></a>
## 2.  Download and set-up

Create a new Python3 environment and install dependencies using `venv` and `pip` or `conda` in terminal:

```
# create and activate an environment using venv
python3 -m venv EpitopeScan_env
source EpitopeScan_env/bin/activate

# or create and activate an environment using conda
conda create --name EpitopeScan_env
conda activate EpitopeScan_env

# then download required packages using pip
pip install -r requirements.txt

# or conda
conda install -c bioconda --file requirements.txt
```
More tutorials on [pip and virtual environments](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) and [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

Clone EpitopeScan repository to your local system by executing: 

```
git clone git@github.com:Aleksandr-biochem/EpitopeScan.git
```

Before using the package, you may run test script as described in [section 3.4](#sec34) to check the performance of EpitopeScan command-line tool.

If you want to call EpitopeScan quickly from anywhere on the system, consult these tutorials on how to add EpitopeScan folder to your system paths on [Windows](https://correlated.kayako.com/article/40-running-python-scripts-from-anywhere-under-windows), [Mac](https://stackoverflow.com/questions/22465332/setting-path-environment-variable-in-osx-permanently) or [Linux](https://www.geeksforgeeks.org/run-python-script-from-anywhere-in-linux/). 

<a name="sec3"></a>
## 3. EpitopeScan usage

<a name="sec31"></a>
### 3.1 Input data

**EpitopeScan requires the following inputs:**

**a)** ***Peptide*** for analysis can be provided as a sequence OR as parent protein name + start and end residue indices. Multiple sequences can be provided as a file (see usage examples below for more).

**b)** ***Genome Multiple Sequence Alignment (MSA)*** in FASTA format. EpitopeScan was originally developed to analyse genome alignments from [COG-UK](https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/). The tool itself does not perform genome alignments. If you want to prepare your set of genomes for analysis, you shoulf use `EPI_ISL_402124.fasta` as the reference sequence and any aligner of your choice (for example, MAFFT as described [here](https://github.com/nihcompmed/SARS-CoV-2-genome)).

**c)** ***Metadata*** on samples dates and lineages in CSV format. Analysis can be performed without this input. However, it will only be possible to count mutations without any insights from sampling date and lineage. EpitopeScan was originally configured to deal with metadata from [COG-UK](https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/). If you wish to use your custom table, make sure to prepare CSV table with the following columns: sequence\_name, sample\_date, epi\_week, usher\_lineage (view `test/example_data/example_metadata.csv` as a guide). 

<a name="sec32"></a>
### 3.2 Running command line tool

`EpitopeScan.py` operates in two modes:
- *scan* to perform MSA file analysis and generate mutation data
- *stat* to generate brief mutation statistics from preexisting mutation data 

The following execution examples are given under the assumption of running the tools from `EpitopeScan` folder as `./EpitopeScan_tool_name.py`. In general, the tools should be called as `path/to/EpitopeScan_tool_name.py`.

<a name="sec321"></a>
#### 3.2.1 Scan mode

*Scan* mode accepts epitope(s), genome MSA and metadata files to generate mutation data. Terminal log messages track run configurations and process (for example, how epitope(s) are mapped onto reference genome and how many genomes are processed at the moment). 

In the end of the run, mutation summary is printed to terminal stdout. It includes the number of mutations for each analysed epitope and the list of discovered mutations with their counts and BLOSUM scores.

Access help section to navigate flags:

```
./EpitopeScan.py scan -h
```

**Required flags:**
- `-e` (--epitope) Single input peptide epitope. Can be specified as name and sequence, comma-separated ("S1,VLLPL"). Or as a peptide name, name of the parent protein with the indices of first and last residue ("S1,S,6,10" is equal to the previous input, protein indexing starts with 1). SARS-CoV-2 protein names can be looked up in `reference_sequences/protein_sequences_reference.fasta`.
- `-f` (--file) A path to file with multiple input peptides. Peptide sequences should be provided in FASTA format. File can also include alternative inputs via residue indices. For example:
```
>S1
VLLPL
>S2,NSP12,7,20
>S3
DYKHYTPSFK
```
- `--msa` Path to input MSA fasta file

**Optional flags:**
- `--metadata` Metadata csv file to merge with mutaion data (absence of metadata limits insights from mutation data)
- `-o` (--out) Output directory name. Optional, default directory name with timestamp is generated automatically
- `-t` (--tag) Sample name pattern to filter. This input string will be compiled as Python regex
- `-q` (--quality\_filter) Upper threshold of max N bases proportion in genome. This proportion is calculated as: count('N' bases in genome without '-' symbols)/length(genome without '-' symbols). Recommended value 0.05
- `-n` (--ambiguity\_threshold) Maximum proportion of ambiguous residues in sample peptide sequence, which is regarded as sufficient coverage. Default 1/3. If  (count(ambiguous bases in peptide)/length(peptide)) > threshold, then sample reported as insufficient coverage.
- `-b` (--blosum) BLOSUM matrix version for mutation scoring. Default 90. See [blosum python package documentation](https://github.com/not-a-feature/blosum) for available options.
- `-s` (--sort), options: {0,1}, Sort mutations summary table by count(0) or score(1). Default 0
- `-a` (--stat), options: {0,1}, Stat individual mutations (0) or combinations(1). Default 0
- `--stat_with_metadata` Only stat samples with metadata in final summary

Basic *scan* run to generate mutation data for 1 peptide and save output to an automatically named folder:

```
./EpitopeScan.py scan -e S1,LTGIAVEQDK --msa test/example_data/example_genomes.fasta
```

Perform analysis with peptide file input and combine mutation data with samples metadata, save to a new folder with custom name:

```
./EpitopeScan.py scan -f test/example_data/example_epitope.fasta --msa test/example_data/example_genomes.fasta --metadata test/example_data/example_metadata.csv -o My_EpitopeScan_Output
```

Add genome quality threshold of 0.07 (7%) and mark any sample with ambiguous bases in epitope's region as insufficient coverage sample:

```
./EpitopeScan.py scan -f test/example_data/example_epitope.fasta --msa test/example_data/example_genomes.fasta --metadata test/example_data/example_metadata.csv -q 0.07 -n 0.0
```

Only keep samples with "England" or "Scotland" mentioned in sequence name:
```
./EpitopeScan.py scan -f test/example_data/example_epitope.fasta --msa test/example_data/example_genomes.fasta --metadata test/example_data/example_metadata.csv -t "England|Scotland"
```

<a name="sec322"></a>
#### 3.2.2 Stat mode

The *stat* mode accepts preexisting *scan* output and prints summary report to terminal stdout with specified options.
There is a help section to navigate flags:

```
./EpitopeScan.py stat -h
```

**Flag description:**
- `-i` (--input) PAth to directory with EpitopeScan scan output
- `-b` (--blosum) BLOSUM version for mutation scoring. Default 90
- `-s` (--sort), options: {0,1}, Sort mutations summary by count(0) or score(1). Default 0
- `-a` (--stat), options: {0,1}, Stat individual mutations (0) or combinations(1).  Default 0
- `--stat_with_metadata`  Only stat samples with metadata
- `--start_date` Subset after this date, dd/mm/yyyy
- `--end_date` Subset before this date, dd/mm/yyyy

Basic run can be launched as follows (you can generate an output from `example_data` to try):

```
./EpitopeScan.py stat -i path/to/EpitopeScan_output_dir
```

To stat occurring combinations instead of individual mutations in a desired date range, and sort summary table by score:

```
./EpitopeScan.py stat -i path/to/EpitopeScan_output_dir -s 1 -a 1 --start_date 15/01/2020 --end_date 17/04/2020
```
<a name="sec323"></a>
#### 3.2.3 Mutation data output

For each peptide, 3 tables are generated and named under following templates:
- `EpitopeName_mutation_data.tsv` contains following columns: sequence\_name, NA\_mutations (list of nucleic acid substitutions, comma-separated), AA\_mutations (list of amino acid substitutions, comma-separated), one-hot-encoding columns detailing presence of discovered AA mutations in each sample, sample\_date, epi\_week (epidemic week), usher\_lineage (viral lineage), has\_metadata (0 if False and 1 if True). If no metadata was provided, the last columns will be filled with NaNs.
- `EpitopeName_AA_mutation_matrix.csv` count matrix. Index corresponds to AA residue in reference epitope sequence. Columns contain counts of each possible residue (including translation stop) and deletions (Δ) at corresponding peptide position
- `EpitopeName_NA_mutation_matrix.csv` count matrix for coding nucleic acid sequece is orginised in same manner to AA matrix

When multiple peptides are provided, the output folder will contain separate subfolders with tables for each peptide.

**Note:** if renamed, the output files will not be suitable for analysis with `EputopeScan stat` or `EpitopeScanGUI`.

<a name="sec33"></a>
### 3.3 Running graphical user-interface tool

EpitopeScan GUI runs locally using the default web-browser as an interface. To launch the application execute:

```
streamlit run ./EpitopeScanGUI.py --server.maxUploadSize 1000
```

Upload mutation data files (one peptide at a time) and explore the interactive summary.
Note that by default streamlit limits file upload to 200 MB. You might want to launch the app specifying a larger limit through `--server.maxUploadSize` depending on the size of your output tables.

<a name="sec34"></a>
### 3.4 Running test scripts

Folder `test` contains `run_EpitopeScan_test.py` script and test data subfolders. This script will perform the analysis for test epitopes and compare EpitopeScan output to the reference mutation data.
Run the test script after installing EpitopeScan to verify the correct performance of the tool:

```
./test/run_EpitopeScan_test.py
```
<a name="sec4"></a>
## 4. Algorithm details

This section contains several remarks about EpitopeScan algorithm, which may be important to correctly interpret the output:

