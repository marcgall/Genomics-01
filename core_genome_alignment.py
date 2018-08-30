#!/usr/bin/python
#
# Take core genome in the fasta format, listing genes by gene family, and
# produces clustalO alignment.
#
# Marc Lorentzen, November 2017.


from Bio import SeqIO
import re
import os
import subprocess
import argparse
from time import gmtime, strftime

parser = argparse.ArgumentParser(description="Takes a core genome from MaGe, aligns gene by gene with clutalO and produces fully aligned core genomes")
parser.add_argument("-i", "--input", metavar="", required=True,
    help="Input core genome, .fasta format. The order of gene families cannot be mixed. Organism strain names must be with spaces as in: '[Oenococcus oeni strain_XXX]'")
parser.add_argument("-o", "--output", metavar="", required=False,
    help="Output file, .fasta format.")
parser.add_argument("-c", "--clustal", metavar="", required=False,
    help="List of arguments to pass to ClustalO. (Remember quotes around the command).")
parser.add_argument("-cl", "--clean", metavar="", type = int, required=False,
    help="Set to 1 to automatically remove all temp files before clustal alignment.\n set to 2 to remove all temporary files. (default behavior)")
parser.add_argument("-x", "--excluded_strains", metavar="", required=False,
    help="A list of strains to exclude in the alignment, separated by whitespace")

args = parser.parse_args()

def get_family_ID(fasta_desc):
    family_ID_temp = re.findall("^\d*\|", fasta_desc)
    family_ID = family_ID_temp[0][:-1]
    return family_ID

def get_strain_name(fasta_desc):
    strain_name_temp1 = re.findall("\[[A-Za-z]* [A-Za-z]* [A-Za-z0-9_-]*]$",
    fasta_desc)
    strain_name_temp2 = re.findall(" [A-Za-z0-9_-]*]" ,strain_name_temp1[0])
    strain_name = strain_name_temp2[0][1:-1]
    return strain_name

def get_stripped_family_ID(fasta_desc):
    family_ID_temp = re.findall("Gene family \d*\|", fasta_desc)
    family_ID = family_ID_temp[0][12:-1]
    return family_ID

def get_stripped_strain_name(fasta_desc):
    strain_name_temp1 = re.findall("\[[A-Za-z0-9_-]*]$",
    fasta_desc)
    strain_name = strain_name_temp1[0][1:-1]
    return strain_name

def validate_input(input_file):
    """
    Verify that every gene family in the input file has exactly one core gene per strain.
    Build a list of exceptions to be skipped.
    Save list of validations/IDs in log.
    """
    total_ID_list_strains = []
    family_ID_list = []
    strain_names = []
    for seq_record in SeqIO.parse(input_file, "fasta"):
        family_ID = get_family_ID(seq_record.description)
        strain_name = get_strain_name(seq_record.description)
        if family_ID not in family_ID_list:
            family_ID_list.append(family_ID)
        if strain_name not in strain_names:
            strain_names.append(strain_name)
        total_ID_list_strains.append((family_ID, strain_name))

    #Find duplicates:
    seen = []
    duplicates = []
    for pair in total_ID_list_strains:
        if pair not in seen:
            seen.append(pair)
        else:
            duplicates.append(pair)
    skip_family_ID = []
    for family_ID, strain_name in duplicates:
        if family_ID not in skip_family_ID:
            skip_family_ID.append(family_ID)
    if skip_family_ID:
        print("More than one gene per strain detected in following gene families (excluded from alignment):")
        for family in skip_family_ID:
            print "Family ID:", family
    else:
        print "No duplicate entries detected."
    #Check that every family has 1 corresponding hit per strain.
    #This is probably a redundant check - but I'll happily sacrifice a minute of calc to be sure.
    missing_family_strain = []
    missing_family = []
    for family in family_ID_list:
        for strain in strain_names:
            if (family, strain) not in total_ID_list_strains:
                missing_family_strain.append((family, strain))
    if missing_family_strain:
        print "Missing gene for strain:"
        for family, strain in missing_family_strain:
            print "Family ID:", family, "strain:", strain
            if family not in missing_family:
                missing_family.append(family)
    else:
        print "No missing genes detected."
    skip_family_ID_total = skip_family_ID + missing_family
    with open("log.txt", "a") as f:
        f.write("\nBefore validation:\nGene Families: " + str(len(family_ID_list)) +
            " Strains: " + str(len(strain_names)) +
            "\nFamilies with more than one entry per strain:\n" + str(skip_family_ID) +
            "\nFamilies with less than one entry per strain:\n" + str(missing_family))
    return skip_family_ID_total

def check_file_add_title(strain_name = "Oenococcus_oeni_xxxx",
    file_name = "Oenococcus_oeni_xxxx.fasta"):
    """
    Check if file is present. If not, create a new one and add the strain name
    on first line.
    """
    if os.path.isfile(file_name):
        pass
    else:
        with open(file_name, "w") as f:
            f.write(">" + strain_name + "\n")

def strip_strain_name(input_string = ">[Oenococcus oeni XXX]"):
    """
    Use the hard-coded Oenococcus regex to strip the strain name from all
    other text.
    """
    strain_name = re.findall("\[Oenococcus.*]", seq_record.description)
    stripped_name = strain_name[0][17:-1]
    stripped_name_underscore = stripped_name.replace(" ", "_")
    return stripped_name, stripped_name_underscore

def get_family_strain_lists(source_fasta, skip_family_ID):
    """
    Iterates through the input fasta file, counts the number of gene objects and
    makes a list of all unique strain names.
    If an object is part of a family excluded in validation, it is skipped.
    """
    #gene_object_count = 0
    unique_family_ID_list = []
    unique_strain_list = []
    for seq_record in SeqIO.parse(source_fasta, "fasta"):
        family_ID = get_family_ID(seq_record.description)
        if family_ID in skip_family_ID:
            continue #skipping to next seq_record
        if family_ID not in unique_family_ID_list:
            unique_family_ID_list.append(family_ID)
        strain_name = get_strain_name(seq_record.description)
        if strain_name not in unique_strain_list:
            unique_strain_list.append(strain_name)
    return unique_family_ID_list, unique_strain_list

def clean_temp_files(clean, unique_strain_list):
    """
    Set to remove the intermediary files between the source fasta file and the
    resulting clustal alignment.
    Note: 'rm *' is not used because the amount of files can become too great
    to handle in one argument.
    """
    if clean > 0:
        print "Cleaning up temporary files..."
        for i in range(10): #This is a hotfix to avoid too many hits on part*
            command_string = "rm core_genomes.part" + str(i) + "*"
            subprocess.call(command_string, shell=True)
        for strain in unique_strain_list:
            command_string = "rm " + strain + ".core.part*"
            subprocess.call(command_string, shell=True)
            if clean > 1:
                command_string = "rm " + strain + ".core.clustal.fasta"
                subprocess.call(command_string, shell=True)

def sort_strain_genes(source_fasta, unique_strain_list, gene_object_count,
    segments, skip_family_ID, excluded_strains):
    """
    Iterates through the source .fasta and outputs a file for each strain, with
    all corresponding genes. To ease computation, the output files are split
    up into smaller segments.
    """
    current_segment = 1
    counter = 0
    print "Iterating through", source_fasta, "Segments:", segments
    for seq_record in SeqIO.parse(source_fasta, "fasta"):
        family_ID = get_family_ID(seq_record.description)
        if family_ID in skip_family_ID:
            continue #skipping to next seq_record
        strain_name = get_strain_name(seq_record.description)
        if excluded_strains:
            if strain_name in excluded_strains:
                continue #Skipping excluded strains.
        counter += 1
        #The conditions for splitting.
        #(= 1 means that we're in the first number in new gene block).
        #First statement ensures no splitting in the middle of core gene blocks
        #Second statement detects if we are going into the next segment
        if (counter % len(unique_strain_list) == 1 and
            counter >= current_segment * gene_object_count / segments):
            current_segment += 1
        #prep strain name.
        strain_name_underscore = strain_name.replace(" ", "_")
        current_file = (strain_name_underscore + ".core.part" +
            str(current_segment) + ".fasta")
        #Create new file. (The tag ASSUMES that segments = max.)
        fasta_desc = "Gene family " + family_ID + "|[" + strain_name + "]"
        #Maybe add handle for later regex?
        check_file_add_title(fasta_desc, current_file)
        #write sequence to file
        with open(current_file, "a") as f:
            f.write(str(seq_record.seq))
        #add a linebreak at the end of all files.
    for strain in unique_strain_list:
        for part in xrange(segments):
            with open(strain + ".core.part" + str(part+1) + ".fasta", "a") as f:
                f.write("\n")
    print "Iteration complete."

def sort_strain_genes_from_clustal(unique_strain_list, segments = 1):
    """
    Iterates through segmented clustal alignments and outputs one file per
    strain. Note: The name of the input files are currently hardcoded.
    """
    for part in xrange(segments):
        for seq_record in SeqIO.parse("core_genomes.part" + str(part+1) +
            ".clustal.fasta", "fasta"):
            strain_name = get_stripped_strain_name(seq_record.description)
            strain_name_underscore = strain_name.replace(" ", "_")
            current_file = strain_name_underscore + ".core.clustal.fasta"
            #Create file
            family_ID = get_stripped_family_ID(seq_record.description)
            check_file_add_title(strain_name, current_file)
            #write sequence to file
            with open(current_file, "a") as f:
                f.write(str(seq_record.seq))
    #Add a newline to end of all files in prep for concaternation.
    for strain in unique_strain_list:
        with open(strain + ".core.clustal.fasta", "a") as f:
            f.write("\n")

#Removed Segments argument from main script, setting it by default to max.
def main_script(source_core_fasta = "test_set.fa", output_file = "alignments.fasta",
                clustal_args = " -v --threads=8", clean = 0, excluded_strains = []):
    """
    Validates core genome calculation output from MaGe (in fasta format).
    Splits the genes into separate files corresponding to each strain,
    then concaternates into a single file, or several segments to ease
    computation.
    Aligns segment of genes in ClustalO,
    """
    print("Run started at " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) +
        "\nValidating input file...")
    with open("log.txt", "w") as f:
        f.write("Log file of core_genome_concat_v7.py run, started " +
            strftime("%Y-%m-%d %H:%M:%S", gmtime()) + "\nInput: " +
            source_core_fasta + "\nOutput: " + output_file + "\nClustal Args: " +
            clustal_args + "\nClean: " + str(clean) + "\nExcluded strains:" + str(excluded_strains))
    skip_family_ID = validate_input(source_core_fasta)
    unique_family_ID_list, unique_strain_list = get_family_strain_lists(source_core_fasta,
        skip_family_ID)
    gene_object_count = len(unique_family_ID_list)*len(unique_strain_list)
    with open("log.txt", "a") as f:
        f.write("\nValidated for run:\n    Number of strains: " + str(len(unique_strain_list)) +
            "\n    Number of core genes: " + str(len(unique_family_ID_list)))
    print("Number of strains: " + str(len(unique_strain_list)) +
        "\nNumber of core genes: " + str(len(unique_family_ID_list)))
    #if segments > len(unique_family_ID_list): #core_gene_count
    segments = len(unique_family_ID_list)
    sort_strain_genes(source_core_fasta, unique_strain_list, gene_object_count,
        segments, skip_family_ID, excluded_strains)
    # Concaternate the groups of strains.
    for part in xrange(segments):
        command_string = ("cat *.core.part" + str(part+1) +
        ".fasta > core_genomes.part" + str(part+1) + ".fasta")
        subprocess.call(command_string, shell=True)
    #Send result to clustalO
        print("Aligning core_genomes.part" + str(part+1) +
            ".fasta with ClustalO.")
        command_string = ("clustalo -i core_genomes.part" + str(part+1) +
            ".fasta -o core_genomes.part" + str(part+1) + ".clustal.fasta " +
            clustal_args)
        subprocess.call(command_string, shell=True)
    print "Building new strain fastas from clustal alignment."
    sort_strain_genes_from_clustal(unique_strain_list, segments)
    #Concaternate aligned strain files.
    command_string = "cat *.core.clustal.fasta > " + output_file
    subprocess.call(command_string, shell=True)
    clean_temp_files(clean, unique_strain_list)
    with open("log.txt", "a") as f:
        f.write("\nRun ended at " + strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    print("Done.")

if __name__ == "__main__":
    #Handling empty args:
    if args.output is None:
        output_file = "core_genome_alignment.fasta"
    else:
        output_file = args.output
    if args.clustal is None:
        clustal_args = "--use-kimura" # --use-kimura
    else:
        clustal_args = args.clustal
    if args.clean is None:
        clean = 2
    else:
        clean = args.clean
    if args.excluded_strains is None:
        excluded_strains = []
    else:
        excluded_strains = args.excluded_strains.split() #Can easily put in different separator
    main_script(args.input, output_file, clustal_args, clean, excluded_strains)
