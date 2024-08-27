#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

fasta_infile = None
min_seq_length = None
output_dir = None

parser.add_argument('--fasta_infile', action='store', dest='fasta_infile',
                    help='input fasta file path as input. (i.e. sequences.fasta)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output fasta file path as input. (i.e. new_sequences.fasta)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_infile = results.fasta_infile
output_dir = results.output_dir

if(fasta_infile == None):
	print('\n')
	print('error: please use the --fasta_infile option to specify the input fasta file path as input')
	print('fasta_infile =' + ' ' + str(fasta_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output dirctory path as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

min_seq_length_ref_db = 1200
max_seq_length_ref_db = 2000

min_seq_length_vsearch_db = 400

basename = os.path.basename(fasta_infile)
filename = os.path.splitext(basename)[0]

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

fasta_outfile1 = os.path.join(output_dir, filename + ".emirge.ref.fasta")
fasta_output_file1 = open(fasta_outfile1, "w+")

fasta_outfile2 = os.path.join(output_dir, filename + ".vsearch.fasta")
fasta_output_file2 = open(fasta_outfile2, "w+")
for record in SeqIO.parse(fasta_infile, "fasta"):
    seq_id = record.id
    desc = record.description
    seq = record.seq
    seq_length = len(seq)
    
    if((int(seq_length) >= int(min_seq_length_ref_db)) and (int(seq_length) <= int(max_seq_length_ref_db))):
        SeqIO.write(record, fasta_output_file1, "fasta")

    if(int(seq_length) >= int(min_seq_length_vsearch_db)):
        SeqIO.write(record, fasta_output_file2, "fasta")

        
fasta_output_file1.close()
fasta_output_file2.close()


