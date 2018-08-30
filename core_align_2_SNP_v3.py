#!/usr/bin/python
#
# Take alignment file. Pick first sequence as the reference to compare all other
# sequences against. Makes a filter sequence of 0/1s based on the logic rules of
# picking: Only non-conserved bases, no N/- characters. The filter is created
# through iteration and is used at the end to filter all the sequences into a
# new output.
#
# Version 2: Output also: A list of the position of the SNPs (in the
# core genome input file).
#
# Version 3: Accepting '-' and N to be used in the output. This means that any
# position with even one missing space will be saved to the output (which may
# not be the best method of finding deletions/insertions).
#
# Marc Lorentzen, November 2017

from Bio import SeqIO
import itertools
import argparse

parser = argparse.ArgumentParser(description="Takes fasta alignment file and removes all conserved bases and N/-'s.")
parser.add_argument("-i", "--input", metavar="", required=True,
    help="Input alignment in fasta format.")
parser.add_argument("-o", "--output", metavar="", required=True,
    help="Output file, .fasta format.")

args = parser.parse_args()

def get_ref_and_filter(input_alignment):
    """
    Get reference strain and initialize the filter sequence.
    (In this version, the first sequence is taken as reference.)
    """
    #Get reference strain:
    ref_seq = []
    for seq_record in SeqIO.parse(input_alignment, "fasta"):
        ref_seq = list(seq_record.seq)
        break
    #Creating the initial state of the filter.
    filter_seq = [0 for i in xrange(len(ref_seq))]
    return ref_seq, filter_seq

def compare_seqs(ref_seq, query_seq, filter_seq):
    """
    Compare base by base of ref and query. Rule for filtering: Remove if base position is conserved
    OR if position is missing '-'.
    """
    new_filter_seq = []
    for ref, query, filt in itertools.izip(ref_seq, query_seq, filter_seq):
        # First find any unwanted characters. Then scan through to find not-conserved positions.
        #if ref in ("N", "-") or query in ("N", "-"):
        #    new_filter_seq.append(2) #2 is here a stand-in to be stripped at end.
        if filt == 0 and ref != query :
            new_filter_seq.append(1)
        else:
            new_filter_seq.append(int(filt))
    return new_filter_seq

def filter_query(query_seq, filter_seq):
    """
    Filters ref_seq using the filter_seq.
    Filter must be list of integers/booleans, not string.
    """
    filtered_seq = list(itertools.compress(query_seq, filter_seq))
    return filtered_seq

def iterate_seqs(input_alignment, output_file):
    """
    The main script. Takes input alignment, gets reference and filter. Iterates through input file to update filter,
    then uses the updated filter on each sequence in turn to produce the filtered alignment output.
    """
    ref_seq, filter_seq = get_ref_and_filter(input_alignment)
    #Iterate through the sequences, updating the filter.
    for seq_record in SeqIO.parse(input_alignment, "fasta"):
        filter_seq = compare_seqs(ref_seq, seq_record.seq, filter_seq)
    #Setting all the '2' elements to 0.
    #filter_seq = [0 if elem == 2 else elem for elem in filter_seq]
    #Use the filter to generate a new file.
    for seq_record in SeqIO.parse(input_alignment, "fasta"):
        filtered_seq = "".join(filter_query(seq_record.seq, filter_seq))
        with open(output_file, "a") as f:
            f.write(">" + seq_record.description + "\n" + filtered_seq + "\n")
    #Get list of SNP positions.
    pos_counter = 0
    pos_list = []
    for pos in filter_seq:
        if pos:
            pos_list.append(pos_counter)
        pos_counter += 1
    with open(output_file + ".poslist", "a") as f:
        for pos in pos_list:
            f.write((str(pos) + "\n"))

if __name__ == "__main__":
    iterate_seqs(args.input, args.output)
    pass
