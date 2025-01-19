# Get the source ~/.bashrc file so that we can use conda.
source ~/.bashrc

# Filtering chimeras.
conda activate vsearch_env

emirge_rename_file="/home/AGR.GC.CA/muirheadk/hybridization_analysis/16S_analysis_output/emirge/16SBC01_S1/16SBC01_S1_emirge.fasta"

# Emirge silva fasta database.
emirge_fasta_db="/home/AGR.GC.CA/muirheadk/hybridization_analysis/databases/vsearch_gtdb/bac120_ssu_reps_r220.fasta"

# The number of cpu threads to use.
num_threads=40

fastq_filename="16SBC01_S1"

emirge_sample_output_dir="/home/AGR.GC.CA/muirheadk/hybridization_analysis/test"
mkdir -p ${emirge_sample_output_dir}

taxonomy_best_hits_file="${emirge_sample_output_dir}/${fastq_filename}_vsearch_taxonomy_best_hits.tsv"
taxonomy_raw_hits_file="${emirge_sample_output_dir}/${fastq_filename}_vsearch_taxonomy_raw.tsv"

vsearch_tax_log_file="${emirge_sample_output_dir}/${fastq_filename}_vsearch_tax_log.txt"

echo "vsearch --usearch_global ${emirge_renamed_file} --db ${emirge_fasta_db} --notrunclabels --userout ${taxonomy_best_hits_file} --userfields query+target+id --uc ${taxonomy_raw_hits_file} --id 0.97 --iddef 0 --log ${vsearch_tax_log_file} --threads ${num_threads} --uc_allhits --maxaccepts 30 --top_hits_only --strand both --gapopen '*'"
vsearch --usearch_global ${emirge_renamed_file} --db ${emirge_fasta_db} --notrunclabels --userout ${taxonomy_best_hits_file} --userfields query+target+id --uc ${taxonomy_raw_hits_file} --id 0.97 --iddef 0 --log ${vsearch_tax_log_file} --threads ${num_threads} --uc_allhits --maxaccepts 30 --top_hits_only --strand both --gapopen '*'


