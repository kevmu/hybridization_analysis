#!/bin/bash
#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G
#SBATCH --output=run_hybrid_emirge_analysis.%A.out
#SBATCH --error=run_hybrid_emirge_analysis.%A.err

# Example to get fastq_files_list.txt

# cpn60
#find /archive/dumonceauxt/230801_M01666_0203_000000000-KTT65_IMPACTTlilacBBSPTWMQhybsEPL/Fastq -type f -name "*.fastq.gz" | grep "cpn60" | rev | cut -d '/' -f1 | rev | sed 's/_L001_R[1-2]_001.fastq.gz//g' | sort -V | uniq > cpn60_fastq_files_list.txt

# Get the source ~/.bashrc file so that we can use conda.
source ~/.bashrc

# The list of fastq filenames.
fastq_list_file="/home/AGR.GC.CA/muirheadk/hybridization_analysis/cpn60_fastq_files_list.txt"

# The fastq input file directory.
fastq_input_dir="/archive/dumonceauxt/230801_M01666_0203_000000000-KTT65_IMPACTTlilacBBSPTWMQhybsEPL/Fastq"
#fastq_input_dir="/home/AGR.GC.CA/muirheadk/hybridization_datasets/cpn60"

# The read1 fastq suffix.
read1_suffix="_L001_R1_001.fastq.gz"

# The read2 fastq suffix.
read2_suffix="_L001_R2_001.fastq.gz"

# The number of cpu threads to use.
num_threads=40

## Emirge databases

# cpn60 UT databases.

# Emirge cpn60 UT fasta database.
emirge_fasta_db="/home/AGR.GC.CA/muirheadk/hybridization_analysis/databases/emirge_dbs/cpn60_db/cpn60_all_nut_seq.clustered.emirge.ref.fasta"

# Emirge cpn60 UT bowtie database.
emirge_bowtie_db_index="/home/AGR.GC.CA/muirheadk/hybridization_analysis/databases/emirge_dbs/cpn60_db/cpn60_all_nut_seq.clustered.emirge.ref"

## Bowtie2 database.

# Bowtie2 database index for mapping reads and automatically determining the parameters for emirge.

# The bowtie2 database index.

bowtie2_db_index="/home/AGR.GC.CA/muirheadk/hybridization_analysis/databases/bowtie2_dbs/cpn60_db/cpn60_all_nut_seq.clustered.emirge.ref"

# The base output directory.
output_dir="/home/AGR.GC.CA/muirheadk/hybridization_analysis/cpn60_analysis_output"
mkdir -p $output_dir

# The preprocessing output directory.
preprocessing_output_dir="${output_dir}/preprocessing"
mkdir -p $preprocessing_output_dir

# The emirge output directory.
emirge_output_dir="${output_dir}/emirge"

### Cutadapt program parameters.

# The hybridization probe primers for cutadapt.

# H279 5'-GAIIIIGCIGGIGAYGGIACIACIAC-3'	
fwd_adapter_H279="GANNNNGCNGGNGAYGGNACNACNAC"

# H280 5'-YKIYKITCICCRAAICCIGGIGCYTT-3'
rev_adapter_H280="YKNYKNTCNCCRAANCCNGGNGCYTT"

# H1612 5'-GAIIIIGCIGGYGACGGYACSACSAC-3'
fwd_adapter_H1612="GANNNNGCNGGYGACGGYACSACSAC"

# H1613 5'-CGRCGRTCRCCGAAGCCSGGIGCCTT-3'
rev_adapter_H1613="CGRCGRTCRCCGAAGCCSGGNGCCTT"


### Prinseq program parameters.

min_qual_mean=25

trim_qual_left=20

trim_qual_right=20

trim_qual_window=3

trim_qual_step=1

trim_qual_type="mean"

trim_qual_rule="lt"

lc_method="dust"

lc_threshold=7

### Cutadapt and Prinseq shared parameters.

# The maximum number of Ns allowed in the reads for filtering.
ns_max_p=15

# The minimum length of fastq reads after filtering.
min_len=60

### EMIRGE parameters.

# ITERATIONS
# Number of iterations to perform.  It may be necessary
# to use more iterations for more complex samples
# (default=40)
num_iter=40

# SNP_FRACTION_THRESH
# If fraction of variants in a candidate sequence
# exceeds this threhold, then split the candidate into
# two sequences for next iteration.  See also
# --variant_fraction_thresh. (default: 0.04)
#snp_fraction_thresh=0.30
snp_fraction_thresh=0.04

# VARIANT_FRACTION_THRESH
# minimum probability of second most probable base at a
# site required in order to call site a variant.  See
# also --snp_fraction_thresh.  (default: 0.10)
#variant_fraction_thresh=0.10
variant_fraction_thresh=0.10

# JOIN_THRESHOLD
# If two candidate sequences share >= this fractional
# identity over their bases with mapped reads, then
# merge the two sequences into one for the next
# iteration.  (default: 0.97; valid range: [0.95, 1.0] )
join_threshold=0.97

# MIN_DEPTH
# minimum average read depth below which a candidate
# sequence is discarded for next iteration (default: 3)
min_depth=3


# Cutadapt


# Activate the cutadapt conda environment.
conda activate cutadapt_env

# Make the cutadapt output directory.
cutadapt_dir="${preprocessing_output_dir}/cutadapt"
mkdir -p $cutadapt_dir

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";

    fastq_read1_file="${fastq_input_dir}/${fastq_filename}${read1_suffix}"
    fastq_read2_file="${fastq_input_dir}/${fastq_filename}${read2_suffix}"
   
    cutadapt_sample_dir="${cutadapt_dir}/${fastq_filename}"
    mkdir -p $cutadapt_sample_dir

    trimmed_fastq_read1_file="${cutadapt_sample_dir}/${fastq_filename}_trimmed_R1.fastq"
    trimmed_fastq_read2_file="${cutadapt_sample_dir}/${fastq_filename}_trimmed_R2.fastq"

    if [ ! -s $trimmed_fastq_read1_file ] && [ ! -s $trimmed_fastq_read2_file ];
    then
        echo "cutadapt -m ${min_len} --max-n ${ns_max_p} -g ${fwd_adapter_H279} -g ${fwd_adapter_H1612} -a ${rev_adapter_H280} -a ${rev_adapter_H1613} -j ${num_threads} -o ${trimmed_fastq_read1_file} -p ${trimmed_fastq_read2_file} ${fastq_read1_file} ${fastq_read2_file}"
        cutadapt -m ${min_len} --max-n ${ns_max_p} -g ${fwd_adapter_H279} -g ${fwd_adapter_H1612} -a ${rev_adapter_H280} -a ${rev_adapter_H1613} -j ${num_threads} -o ${trimmed_fastq_read1_file} -p ${trimmed_fastq_read2_file} ${fastq_read1_file} ${fastq_read2_file}
    else

        trimmed_fastq_read1_filename=$(basename $trimmed_fastq_read1_file)
        trimmed_fastq_read2_filename=$(basename $trimmed_fastq_read2_file)
        echo "The following files have already been created."
        echo "${trimmed_fastq_read1_filename}"
        echo "${trimmed_fastq_read2_filename}"
        echo "Skipping to next set of commands!!!"

    fi

done

### Prinseq 

# Activate the prinseq conda environment.
conda activate prinseq_env

# Make the prinseq output directory.
prinseq_dir="${preprocessing_output_dir}/prinseq"
mkdir -p $prinseq_dir

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";

    cutadapt_sample_dir="${cutadapt_dir}/${fastq_filename}"
    mkdir -p $cutadapt_sample_dir

    trimmed_fastq_read1_file="${cutadapt_sample_dir}/${fastq_filename}_trimmed_R1.fastq"
    trimmed_fastq_read2_file="${cutadapt_sample_dir}/${fastq_filename}_trimmed_R2.fastq"

    prinseq_sample_dir="${prinseq_dir}/${fastq_filename}"
    mkdir -p $prinseq_sample_dir

    fastq_file_prefix="${prinseq_sample_dir}/${fastq_filename}_filtered"

    filtered_fastq_read1_file="${prinseq_sample_dir}/${fastq_filename}_filtered_1.fastq"
    filtered_fastq_read2_file="${prinseq_sample_dir}/${fastq_filename}_filtered_2.fastq"

    if [ ! -s $filtered_fastq_read1_file ] && [ ! -s $filtered_fastq_read2_file ];    
    then
        echo "perl software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ${trimmed_fastq_read1_file} -fastq2 ${trimmed_fastq_read2_file} -out_format 3 -out_good ${fastq_file_prefix} -out_bad null -no_qual_header -min_qual_mean ${min_qual_mean} -ns_max_p ${ns_max_p} -lc_method ${lc_method} -lc_threshold ${lc_threshold} -trim_qual_left ${trim_qual_left} -trim_qual_right ${trim_qual_right} -trim_qual_type ${trim_qual_type} -trim_qual_rule ${trim_qual_rule} -trim_qual_window ${trim_qual_window} -trim_qual_step ${trim_qual_step} -min_len ${min_len} -verbose"
    	perl software/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ${trimmed_fastq_read1_file} -fastq2 ${trimmed_fastq_read2_file} -out_format 3 -out_good ${fastq_file_prefix} -out_bad null -no_qual_header -min_qual_mean ${min_qual_mean} -ns_max_p ${ns_max_p} -lc_method ${lc_method} -lc_threshold ${lc_threshold} -trim_qual_left ${trim_qual_left} -trim_qual_right ${trim_qual_right} -trim_qual_type ${trim_qual_type} -trim_qual_rule ${trim_qual_rule} -trim_qual_window ${trim_qual_window} -trim_qual_step ${trim_qual_step} -min_len ${min_len} -verbose
    else
    	filtered_fastq_read1_filename=$(basename $filtered_fastq_read1_file)
    	filtered_fastq_read2_filename=$(basename $filtered_fastq_read2_file)
	echo "The following files have already been created."
	echo "${filtered_fastq_read1_filename}"
	echo "${filtered_fastq_read2_filename}"
	echo "Skipping to next set of commands!!!"
    fi

done

# Make the insert size output directory.
insert_size_dir="${output_dir}/insert_size"
mkdir -p $insert_size_dir

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";

    prinseq_sample_dir="${prinseq_dir}/${fastq_filename}"
    mkdir -p $prinseq_sample_dir

    filtered_fastq_read1_file="${prinseq_sample_dir}/${fastq_filename}_filtered_1.fastq"
    filtered_fastq_read2_file="${prinseq_sample_dir}/${fastq_filename}_filtered_2.fastq"
        
    insert_size_sam_file="${insert_size_dir}/${fastq_filename}_insert_size.sam"
    insert_size_bam_file="${insert_size_dir}/${fastq_filename}_insert_size.bam"
    insert_size_stats_file="${insert_size_dir}/${fastq_filename}_insert_size_stats.txt"

    # Activate the bowtie2 conda environment. 
    conda activate bowtie2_env
    
    if [ ! -s $insert_size_sam_file ];
    then
        echo "bowtie2 -x ${bowtie2_db_index} -p ${num_threads} -1 ${filtered_fastq_read1_file} -2 ${filtered_fastq_read2_file} > ${insert_size_sam_file}"
        bowtie2 -x ${bowtie2_db_index} -p ${num_threads} -1 ${filtered_fastq_read1_file} -2 ${filtered_fastq_read2_file} > ${insert_size_sam_file}
    else
        insert_size_sam_filename=$(basename $insert_size_sam_file)
        echo "The ${insert_size_sam_filename} file has already been created. Skipping to next set of commands!!!"
    fi
    
    # Activate the samtools conda environment.
    conda activate samtools_env
    
    if [ ! -s $insert_size_bam_file ];
    then
        echo "samtools view -bS < ${insert_size_sam_file} > ${insert_size_bam_file}"
        samtools view -bS < ${insert_size_sam_file} > ${insert_size_bam_file}
    else
        insert_size_bam_filename=$(basename $insert_size_sam_file)
        echo "The ${insert_size_bam_filename} file has already been created. Skipping to next set of commands!!!"
    fi
    
    if [ ! -s $insert_size_stats_file ];
    then
        echo "samtools stats ${insert_size_bam_file} > ${insert_size_stats_file}"
        samtools stats ${insert_size_bam_file} > ${insert_size_stats_file}
    else
        insert_size_stats_filename=$(basename $insert_size_stats_file)
        echo "The ${insert_size_stats_filename} file has already been created. Skipping to next set of commands!!!"
    fi
    
done



### Emirge

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";
    
    prinseq_sample_dir="${prinseq_dir}/${fastq_filename}"
    mkdir -p $prinseq_sample_dir

    filtered_fastq_read1_file="${prinseq_sample_dir}/${fastq_filename}_filtered_1.fastq"
    filtered_fastq_read2_file="${prinseq_sample_dir}/${fastq_filename}_filtered_2.fastq"
 
    insert_size_stats_file="${insert_size_dir}/${fastq_filename}_insert_size_stats.txt"
    
    # MAX_READ_LENGTH
    # length of longest read in input data.
    #SN      maximum length: 151
    max_read_length=$(grep "maximum length" ${insert_size_stats_file} | cut -f3 | xargs printf "%.0f\n")
    
    # INSERT_MEAN
    # insert size distribution mean.
    #SN      insert size average:    227.2
    # Get the rounded value for the insert size average because emirge needs a whole number.
    insert_mean=$(grep "insert size average" ${insert_size_stats_file} | cut -f3 | xargs printf "%.0f\n")
    
    # INSERT_STDDEV
    # insert size distribution standard deviation.
    #SN      insert size standard deviation: 79.0
    # Get the rounded value for the insert size standard deviation because emirge needs a whole number.
    insert_stddev=$(grep "insert size standard deviation" ${insert_size_stats_file} | cut -f3 | xargs printf "%.0f\n")
       

    # Run emirge

    # The emirge sample output directory. 
    emirge_sample_output_dir="${emirge_output_dir}/${fastq_filename}"
    mkdir -p ${emirge_sample_output_dir}

    # Activate the emirge conda environment.
    conda activate emirge_env

    # The emirge fasta file.
    emirge_fasta_file="${emirge_sample_output_dir}/iter.${num_iter}/iter.${num_iter}.cons.fasta"

    # Check if the emirge fasta file is created. If not then run emirge otherwise skip to next set of commands.
    if [ ! -s $emirge_fasta_file ];
    then

        echo "emirge.py ${emirge_sample_output_dir} -1 ${filtered_fastq_read1_file} -2 ${filtered_fastq_read2_file} -f ${emirge_fasta_db} -b ${emirge_bowtie_db_index} -l ${max_read_length} -i ${insert_mean} -s ${insert_stddev} -n ${num_iter} -a ${num_threads} -p  ${snp_fraction_thresh} -v ${variant_fraction_thresh} -j ${join_threshold} -c ${min_depth} --phred33"
        emirge.py ${emirge_sample_output_dir} -1 ${filtered_fastq_read1_file} -2 ${filtered_fastq_read2_file} -f ${emirge_fasta_db} -b ${emirge_bowtie_db_index} -l ${max_read_length} -i ${insert_mean} -s ${insert_stddev} -n ${num_iter} -a ${num_threads} -p  ${snp_fraction_thresh} -v ${variant_fraction_thresh} -j ${join_threshold} -c ${min_depth} --phred33

    else
        emirge_fasta_filename=$(basename $emirge_fasta_file)
        echo "The ${emirge_fasta_filename} file has already been created. Skipping to next set of commands!!!"
    fi


    emirge_iter_dir="${emirge_sample_output_dir}/iter.${num_iter}"
    emirge_renamed_file="${emirge_sample_output_dir}/${fastq_filename}_emirge.fasta"

    if [ ! -s $emirge_renamed_file ];
    then

	echo "emirge_rename_fasta.py ${emirge_iter_dir} > ${emirge_renamed_file}"
    	emirge_rename_fasta.py ${emirge_iter_dir} > ${emirge_renamed_file}
    
    else
        emirge_renamed_filename=$(basename $emirge_renamed_file)
        echo "The ${emirge_renamed_filename} file has already been created. Skipping to next set of commands!!!"
    fi

done

for fastq_filename in $(cat $fastq_list_file);
do
	echo "Processing ${fastq_filename}......";
    	
	emirge_sample_output_dir="${emirge_output_dir}/${fastq_filename}"
        mkdir -p ${emirge_sample_output_dir}

	emirge_renamed_file="${emirge_sample_output_dir}/${fastq_filename}_emirge.fasta"
	taxonomy_best_hits_file="${emirge_sample_output_dir}/${fastq_filename}_vsearch_taxonomy_best_hits.tsv"
	taxonomy_raw_hits_file="${emirge_sample_output_dir}/${fastq_filename}_vsearch_taxonomy_raw.tsv"

	vsearch_tax_log_file="${emirge_sample_output_dir}/${fastq_filename}_vsearch_tax_log.txt"

	conda activate vsearch_env
	
	echo "vsearch --usearch_global ${emirge_renamed_file} --db ${emirge_fasta_db} --notrunclabels --userout ${taxonomy_best_hits_file} --userfields query+target+id --uc ${taxonomy_raw_hits_file} --id 0.97 --iddef 0 --log ${vsearch_tax_log_file} --threads ${num_threads} --uc_allhits --maxaccepts 30 --top_hits_only --strand both --gapopen '*'"
	vsearch --usearch_global ${emirge_renamed_file} --db ${emirge_fasta_db} --notrunclabels --userout ${taxonomy_best_hits_file} --userfields query+target+id --uc ${taxonomy_raw_hits_file} --id 0.97 --iddef 0 --log ${vsearch_tax_log_file} --threads ${num_threads} --uc_allhits --maxaccepts 30 --top_hits_only --strand both --gapopen '*'

done

