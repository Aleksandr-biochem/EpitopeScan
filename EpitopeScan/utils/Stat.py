import os
import glob
import blosum as bl
import pandas as pd

from EpitopeScan.utils.ProteinUtils import *


def StatIndividualMutations(df, blosum_version=90):
    """
    Count and score individual protein mutations

    df - pd DataFrame, mutations data

    Returns:
    AA_mutations_stat - dict(str(mutation): (int(count), float(score)))
    """
    blosum_matrix = bl.BLOSUM(blosum_version)

    # assumes that first 3 columns are sequence name and NA, AA mutation lists
    # and last 4 columns are metadata related
    stat = df.iloc[:, 3:-4].apply(lambda x: (x == 1).sum())

    # create stat dictionary
    AA_mutations_stat = dict()
    for item in stat.items():
        if "Δ" in item[0]:
            score = -6.0  # some negative score for deletion
        else:
            score = blosum_matrix[item[0][0]][item[0][-1]]
        AA_mutations_stat[item[0]] = (item[1], score)

    return AA_mutations_stat


def CombinationsStat(df, blosum_version=90):
    """
    Count and score protein mutations as ocurring combinations

    df - pd DataFrame, mutations data

    Returns:
    AA_mutations_stat - dict(str(mutation): (int(count), float(score)))
    """
    blosum_matrix = bl.BLOSUM(blosum_version)

    # create the stat dictionary
    AA_mutations_stat = dict()

    # get mutation combinations
    all_mutations = set(df.columns[3:-4])
    mutations, combinations = df.iloc[:, 3:-4], []
    if mutations.shape[1] > 0:
        combinations = mutations.drop_duplicates()
        combinations = combinations.apply(lambda row: list(row[row == 1].index), axis=1)

    # count and score combinations
    for combination in combinations:
        if len(combination) > 0:
            mutations_filt = mutations

            for mut in combination:
                mutations_filt = mutations_filt[mutations_filt[mut] == 1]
            for mut in all_mutations - set(combination):
                mutations_filt = mutations_filt[mutations_filt[mut] == 0]

            count = mutations_filt.shape[0]

            score = 0.0
            for mut in combination:
                if "Δ" in mut:
                    score -= 6.0  # subtract deletion score
                else:
                    score += blosum_matrix[mut[0]][mut[-1]]

            AA_mutations_stat[",".join(combination)] = (count, score)

    return AA_mutations_stat


def MutationTextSummary(peptide, df, stat_mutations=0, sort_by=0, blosum_version=90):
    """
    Print protein mutation summary to stdout

    peptide - Protein instance, peptide to summarise
    df - pd DataFrame, mutation data
    stat_mutations - int, 0 to stat individual mutations, 1 - mutation combinations
    sort_by - int, 0 to sort mutation report by count, 1 - by BLOSUM score
    blosum_version - BLOSUM version to use for scoring, default 90

    Returns:
    summary_table - pd DataFrame
    """

    # count and score mutations
    if stat_mutations == 0:
        AA_mutations_stat = StatIndividualMutations(df, blosum_version)
    elif stat_mutations == 1:
        AA_mutations_stat = CombinationsStat(df, blosum_version)

    # sort AA substitutions by counts or scores
    key_to_sort = {
        0: (lambda x: x[1][0]),  # by count
        1: (lambda x: x[1][1]),
    }  # by score
    AA_mutations_stat = dict(
        sorted(AA_mutations_stat.items(), key=key_to_sort[sort_by], reverse=True)
    )

    # remove 0 counts
    AA_mutations_stat = {k: v for k, v in AA_mutations_stat.items() if v[0] > 0}

    # generate summary table
    summary_table = []
    for mut in AA_mutations_stat:
        summary_row = []
        mut_string = ["_"] * len(peptide)

        for m in mut.split(","):
            position = int(m[1:-1])
            AA2 = "-" if ("Δ" in m) else m[-1]
            mut_string[position - peptide.protein_start] = AA2

        summary_row.append("".join(mut_string))
        summary_row.append(f"{mut}")
        summary_row.append(f"{AA_mutations_stat[mut][0]}")
        summary_row.append(f"{AA_mutations_stat[mut][1]}")

        summary_table.append(summary_row)

    # return table as pd DataFrame
    summary_table = pd.DataFrame(
        summary_table, columns=[peptide.sequence, " ", "count", "score"]
    )

    return summary_table


def GetKeyStat(mutation_data, stat_only_metadata):
    """
    Return and print key statistics of loaded samples

    peptide_name - str, peptide name
    mutation_data - pd DataFrame, data to stat
    stat_only_metadata - bool, stat only samples with metadata

    Returns
    key_stat - dict, key samples statistics
    """
    total = mutation_data.shape[0]
    has_metadata = sum(mutation_data["has_metadata"] == 1)

    mask = (
        mutation_data["has_metadata"] == 1
        if stat_only_metadata
        else mutation_data["has_metadata"].isin([0, 1])
    )

    # some_samples_lack_metadata = True if 0 in mutation_data["has_metadata"] else False

    no_mutation = sum((mutation_data["AA_mutations"].isna()) & mask)
    no_coverage = sum((mutation_data["AA_mutations"] == "-") & mask)
    non_functional = sum((mutation_data["AA_mutations"] == "NF") & mask)

    # calculate number of samples with mutations
    have_mutation = has_metadata if stat_only_metadata else total
    have_mutation = have_mutation - no_mutation - no_coverage - non_functional

    # create dictionary
    key_stat = {
        "Have_mutation": have_mutation,
        "No_mutation": no_mutation,
        "No_coverage": no_coverage,
        "Non_functional": non_functional,
        "Total": total,
        "Has_metadata": has_metadata,
    }

    return key_stat


def StatEpitopeData(
    input_dir,
    proteome,
    reference_genome,
    stat_mutations=0,
    sort_by=0,
    time_start=None,
    time_end=None,
    blosum_version=90,
    metadata_filter=False,
):
    """
    Read EpitopeScan output data, stat and print summary

    input_dir - str, path to EpitopeScan results dir
    proteome - list(Protein instances), reference proteome
    reference_genome - str, reference genome sequence
    stat_mitations - int, 0 (stat individual mutations) or 1 (combinations)
    sort_by - int, 0 (by count) or 1 (by BLOSUM score)
    time_start - datetime, start date to stat mutations
    end_time - datetime, end date to stat mutations
    blosum_version - int, BLOSUM version to use for mutation scoring
    metadata_filter - bool, filter samples with metadata
    """

    # check if we work with one or several peptides
    list_dir = [
        os.path.isfile(os.path.join(input_dir, f)) for f in os.listdir(input_dir)
    ]
    one_peptide = all(list_dir)

    # list peptide names to analyse from file names
    if one_peptide:
        # get peptide name from mutation data file
        peptide_names = glob.glob(f"{input_dir}/*.tsv")[0]
        peptide_names = [peptide_names.split("/")[-1].split("_")[0]]
        prefixes = [f"{input_dir}/"]
    else:
        # subdir names are peptide names
        peptide_names = next(os.walk(input_dir))[1]
        prefixes = [f"{input_dir}/{name}/" for name in peptide_names]

    # read data and calculate statistics for each peptide
    printed_headline_statistics, message = False, ""
    for name, prefix in zip(peptide_names, prefixes):
        # extract peptide sequence from AA substitution matrix
        matrix = pd.read_csv(f"{prefix}{name}_AA_mutation_matrix.csv", index_col=0)
        peptide = Protein(name, "".join(matrix.index))

        # map peptide onto proteome
        peptide = MapPeptides([peptide], proteome, reference_genome, verbose=False)[0]

        # read mutation data
        df = pd.read_csv(f"{prefix}{name}_mutation_data.tsv", sep="\t", index_col=0)

        # convert sample_date to datetime
        df["sample_date"] = pd.to_datetime(df["sample_date"], format="%Y-%m-%d")

        # record total number of analysed samples
        if not printed_headline_statistics:
            message += f"Found data on {len(peptide_names)} peptides from {df.shape[0]} samples\n"

        # filter samples with metadata if required
        if metadata_filter:
            df = df[df["has_metadata"] == 1]
            if not printed_headline_statistics:
                message += f"Only keeping {df.shape[0]} samples with metadata\n"

        # if some samples have no metadata
        # date to filters are irrelevent
        some_samples_lack_metadata = (
            False if (sum(df["has_metadata"] == 0) == 0) else True
        )

        # filter by date and get reported range
        if not some_samples_lack_metadata:
            if time_start is not None:
                df = df[df["sample_date"] >= time_start]
                start_date = time_start.strftime("%d/%m/%Y")
            else:
                start_date = df["sample_date"].min().strftime("%d/%m/%Y")

            if time_end is not None:
                df = df[df["sample_date"] <= time_end]
                end_date = time_end.strftime("%d/%m/%Y")
            else:
                end_date = df["sample_date"].max().strftime("%d/%m/%Y")

            message += (
                f"{df.shape[0]} samples dated between {start_date} and {end_date}\n"
            )

        # headline statistics
        if not printed_headline_statistics:
            print(message)
            printed_headline_statistics = True

        # print corresponding summary table
        key_stat = GetKeyStat(df, stat_only_metadata=False)
        print(f"For {peptide.name}")
        total = key_stat["Total"]
        print(
            f"{key_stat['Have_mutation']} ({round(key_stat['Have_mutation'] * 100/ total, 2)}%) have mutations"
        )
        print(
            f"{key_stat['No_mutation']} ({round(key_stat['No_mutation'] * 100/ total, 2)}%) have no mutations"
        )
        print(
            f"{key_stat['No_coverage']} ({round(key_stat['No_coverage'] * 100/ total, 2)}%) reported with insufficient coverage"
        )
        print(
            f"{key_stat['Non_functional']} ({round(key_stat['Non_functional'] * 100/ total, 2)}%) reported as non-functional"
        )
        print()
        summary = MutationTextSummary(
            peptide, df, stat_mutations, sort_by, blosum_version
        )

        print(f"{summary.shape[0]} mutations for {peptide.name}")
        if summary.shape[0] > 0:
            print(summary.to_string())
        print()

    return
