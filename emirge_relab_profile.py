#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from natsort import natsorted

emirge_renamed_fasta_file = "/home/AGR.GC.CA/muirheadk/hybridization_analysis/16S_analysis_output/emirge/16SBC01_S1/16SBC01_S1_emirge.fasta"

output_dir = "/home/AGR.GC.CA/muirheadk/hybridization_analysis/testing_dir"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# The basename of the file.
basename = os.path.basename(emirge_renamed_fasta_file)

# The filename of the file.
filename = os.path.splitext(basename)[0]

csv_writer_file_handle = open(os.path.join(output_dir, "_".join([filename, "profile"]) + ".txt"), "w+")
csv_writer = csv.writer(csv_writer_file_handle, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
csv_writer.writerow(["seq_id","db_id","num_prior","num_norm_prior","seq_length"])
        
emirge_relab_dict = {}            
for record in SeqIO.parse(emirge_renamed_fasta_file, "fasta"):
    seq_id = record.id
    desc = record.description
    seq = record.seq

#    print(desc)

    # 75|CP014031.722667.724218 Prior=0.134884 Length=1523 NormPrior=0.124875
    (seq_id_str, num_prior_str, seq_length_str, num_norm_prior_str) = desc.split(" ")
    print(seq_id_str, num_prior_str, seq_length_str, num_norm_prior_str)

    seq_id = seq_id_str.split("|")[0]
    db_id = seq_id_str.split("|")[1]

    num_prior = num_prior_str.split("=")[1]
    seq_length = seq_length_str.split("=")[1]
    num_norm_prior = num_norm_prior_str.split("=")[1]

    if(not(num_norm_prior in emirge_relab_dict)):
        emirge_relab_dict[str(num_norm_prior)] = []
        emirge_relab_dict[str(num_norm_prior)].append([seq_id, db_id, num_prior, num_norm_prior, seq_length])

    else:
        emirge_relab_dict[str(num_norm_prior)].append([seq_id, db_id, num_prior, num_norm_prior, seq_length])

print(emirge_relab_dict)

# Sort the entries based on the normprior values.
for num_norm_prior in natsorted(emirge_relab_dict, reverse=True):
    
    for row in emirge_relab_dict[num_norm_prior]:
        print(row)

        (seq_id, db_id, num_prior, num_norm_prior, seq_length) = row
        csv_writer.writerow([seq_id, db_id, num_prior, num_norm_prior, seq_length])

csv_writer_file_handle.close()


