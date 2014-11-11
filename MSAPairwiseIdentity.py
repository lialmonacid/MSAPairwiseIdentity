#!/usr/bin/env python

import sys
import getopt
import os
from Bio import SeqIO

OPT_INPUT_FILE = False
OPT_OUTPUT_FILE = False

def Usage():
    print "\nMSAPairwiseIdentity.py is a program that read a multiple sequence alignment in FASTA format and calculate all the pairwise indentities.\n"
    print "Usage:"
    print "\tMSAPairwiseIdentity.py -i [FASTA file]\n"
    print "\nMandatory options:"
    print "\t-i, --input=FILE"
    print "\t\tThe input FASTA file that contains the sequences were the pairwise identity is going to be calculate. "
    print "\nOther options:"
    print "\t-h, --help"
    print "\t\tShow the options of the program."
    print "\t-o, --output=FILE"
    print "\t\tWrite the output to the given file in FASTA format. By default this option is not set and the ouput is written to the STDOUT."
    print "\n"
    sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
    if len(argv) == 0:
        Usage()
    options, remaining = getopt.getopt(argv, 'i:o:h', ['input=','output=','help'])
    opt_flag = {'i': False, 'o':False}
    global OPT_INPUT_FILE, OPT_OUTPUT_FILE
    for opt, argu in options:
        if opt in ('-i', '--input'):
            if not opt_flag['i']:
                if os.path.exists(argu):
                    OPT_INPUT_FILE = argu
                    opt_flag['i'] = True
                else:
                    print >> sys.stderr , "\n[ERROR]: File or path of the input file does not exist. ", argu, "\n"
                    sys.exit(1)
            else:
                print >> sys.stderr , "\n[ERROR]: Trying to redefine the input file. Option -i / --input was already set.\n"
                sys.exit(1)
        elif opt in ('-o', '--output'):
            if not opt_flag['o']:
                if not os.path.dirname(argu): # Empty path means the current directory.
                    OPT_OUTPUT_FILE = argu
                    opt_flag['o'] = True
                else:
                    if os.path.exists(os.path.dirname(argu)):
                        OPT_OUTPUT_FILE = argu
                        opt_flag['o'] = True
                    else:
                        print >> sys.stderr , "\n[ERROR]: Path to write the output does not exist. ", os.path.dirname(argu), "\n"
                        sys.exit(1)
            else:
                print >> sys.stderr , "\n[ERROR]: Trying to redefine the output file. Option -o / --output was already set.\n"
                sys.exit(1)
        elif opt in ('-h', '--help'):
            Usage()
    
    if not opt_flag['i']:
        print >> sys.stderr , "\n[ERROR]: Input file not defined. Option -i / --input.\n"
        sys.exit(1)
    
    if OPT_OUTPUT_FILE: # Setting the output
        OPT_OUTPUT_FILE=open(OPT_OUTPUT_FILE,"w")
    else:
        OPT_OUTPUT_FILE=sys.stdout

def SequenceLengthWithoutGaps(sequence):
    length = 0
    for character in sequence:
        if character != "-":
            length = length + 1
    return length

def CalculateIdenticalAndNonIdenticalCharacters(sequence_1,sequence_2):
    assert len(sequence_1) == len(sequence_2)
    
    count_of_identical_characters = 0
    count_of_non_identical_characters = 0
    count_of_characters_facing_gap = 0
    for i in range(0,len(sequence_1)):
        if sequence_1[i] != "-" and sequence_2[i] != "-" and sequence_1[i] != "X" and sequence_2[i] != "X" and sequence_1[i] == sequence_2[i]:
            count_of_identical_characters = count_of_identical_characters + 1
        elif (sequence_1[i] != "-" and sequence_2[i] != "-" and sequence_1[i] != sequence_2[i]) or (sequence_1[i] == "X" and sequence_2[i] == "X"):
            count_of_non_identical_characters = count_of_non_identical_characters + 1
        elif (sequence_1[i] == "-" and sequence_2[i] != "-") or (sequence_1[i] != "-" and sequence_2[i] == "-"):
            count_of_characters_facing_gap = count_of_characters_facing_gap + 1

    return count_of_identical_characters,count_of_non_identical_characters,count_of_characters_facing_gap

# Parse command line
SetOptions(sys.argv[1:])

msa = []
# Read and load the entire MSA into memory
for record in SeqIO.parse(open(OPT_INPUT_FILE, "rU"), "fasta"):
    msa.append([record.id,SequenceLengthWithoutGaps(record.seq),record.seq])

# Iterate over the MSA and calculate all pairwise identities. A total of ((n*(n-1))/2) comparisons where n is the number of the sequences in the input file.
for i in range(0,len(msa)):
    j = 0
    while j != i:
        number_of_identical_characters,non_identical_characters,characters_facing_gap = CalculateIdenticalAndNonIdenticalCharacters(msa[i][2],msa[j][2])
        pairwise_identity = float(number_of_identical_characters) / float(min(msa[i][1],msa[j][1]))
        
        print >> OPT_OUTPUT_FILE ,msa[i][0] + "\t" + msa[j][0]  + "\t" + "%.8f" %(pairwise_identity) + "\t" + str(non_identical_characters) + "\t" + str(characters_facing_gap)
        j = j + 1

