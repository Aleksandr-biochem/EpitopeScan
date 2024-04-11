#! /usr/bin/env python
import os
import glob
import argparse
from subprocess import Popen

# this script automates Multiple Sequence Alignment update
# with MAFFT

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update MSA with MAFFT aligner")

    parser.add_argument(
        "-s",
        "--sequences_dir",
        help="Path to directory with genome sequences to align",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-a",
        "--alignment_file",
        help="Alignment file to append sequences to",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # list sequence names in alignment
    cur_alignment_file = args.alignment_file
    aligned_sequences = []
    with open(cur_alignment_file, "r") as f:
        aligned_sequences = [l.strip()[1:] for l in f if l.startswith(">")]
    aligned_sequences = set(aligned_sequences)
    print(
        f"Current alignment {args.alignment_file} has {len(aligned_sequences)} sequences"
    )

    # check what sequences are absent from the alignment
    file_stat = {}
    n_seq, n_new_seq = 0, 0
    for file in glob.glob(f"{args.sequences_dir}/*"):
        # check format
        if file.split(".")[-1] not in ["fa", "fasta"]:
            continue

        # read sequence names from file
        with open(file, "r") as f:
            file_sequences = [l.strip()[1:] for l in f if l.startswith(">")]

        # record number of sequences and
        n_seq += len(file_sequences)
        n_new_seq += len(set(file_sequences) - aligned_sequences)
        file_stat[file] = [
            len(file_sequences),
            len(set(file_sequences) - aligned_sequences),
        ]

    print(f"Found {len(file_stat)} files in {args.sequences_dir}")
    print(
        f"Overall {n_seq} sequences with {n_new_seq} sequences absent from current alignment\n"
    )

    # perform alignment of the sequences from new files
    print("Starting the alignment extension...")
    for file in file_stat:
        # if file contains new sequences
        # append them to the alignment
        if file_stat[file][1] > 0:
            print(f"Appending sequences from {file}")

            # if all sequences in a file are new
            if file_stat[file][0] == file_stat[file][1]:
                file_with_sequences = file

            else:
                # create tmp file with new sequences
                with open(f"{file}.tmp", "w") as t:
                    with open(file, "r") as f:
                        new_seq = None
                        for line in f:
                            if line.startswith(">") and (
                                line.strip()[1:] not in aligned_sequences
                            ):
                                new_seq = line
                            elif new_seq is not None:
                                t.write(new_seq)
                                t.write(line)
                                new_seq = None

                file_with_sequences = f"{file}.tmp"

            # run MAFFT
            mafft_args = [
                "mafft",
                "--auto",
                "--keeplength",
                "--addfragments",
                file_with_sequences,
                cur_alignment_file,
            ]
            updated_alignment = open(f"{cur_alignment_file}.upd", "w")
            process = Popen(mafft_args, stdout=updated_alignment)
            returncode = process.wait()
            updated_alignment.close()

            # remove the old alignment file
            os.remove(cur_alignment_file)
            os.rename(f"{cur_alignment_file}.upd", cur_alignment_file)

            # delete tmp file if neccessary
            if ".tmp" in file_with_sequences:
                os.remove(file_with_sequences)
