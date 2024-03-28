#!/bin/bash
#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=100G
#SBATCH --output=run_hybrid_emirge_analysis.%A.out
#SBATCH --error=run_hybrid_emirge_analysis.%A.err


# Get the source ~/.bashrc file so that we can use conda.
source ~/.bashrc

# The fastq input file directory.
fastq_input_dir="/archive/dumonceauxt/230801_M01666_0203_000000000-KTT65_IMPACTTlilacBBSPTWMQhybsEPL/Fastq"

# The read1 fastq suffix.
read1_suffix="_R1.fastq"

# The read2 fastq suffix.
read2_suffix="_R2.fastq"

fastq_file_ext=".fastq"

# The number of cpu threads to use.
num_threads=8

# Emirge databases
emirge_fasta_db="/home/AGR.GC.CA/muirheadk/hybridization_analysis/emirge_db/silva_db/SILVA_138_SSURef_NR99_tax_silva_trunc.fixed.clustered.emirge.ref.fasta"

emirge_bowtie_db="/home/AGR.GC.CA/muirheadk/hybridization_analysis/emirge_db/silva_db/SILVA_138_SSURef_NR99_tax_silva_trunc.fixed.clustered.emirge.ref"

# The base output directory.
output_dir="/home/AGR.GC.CA/muirheadk/hybridization_analysis/impactt_hybrid"
mkdir -p $output_dir

# The preprocessing output directory.
preprocessing_output_dir="${output_dir}/preprocessing"
mkdir -p $preprocessing_output_dir

fastq_list_file="${preprocessing_output_dir}/fastq_files_list.txt"

# The emirge output directory.
emirge_output_dir="${output_dir}/emirge"

### Trimmomatic program parameters.

# Cut bases off the start of a read, if below a threshold quality < N.
cut_leading_bases=3

# Cut bases off the end of a read, if below a threshold quality < N.
cut_trailing_bases=20

## Use a sliding window of size sliding_window_size that will remove bases if their phred score is below sliding_window_phred_score.

# Sliding window of size N.
sliding_window_size=4

# Remove bases if their phred score is below 20
sliding_window_phred_score=15

# Minimum length for fastq reads.
min_read_length=36

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
snp_fraction_thresh=0.30

# VARIANT_FRACTION_THRESH
# minimum probability of second most probable base at a
# site required in order to call site a variant.  See
# also --snp_fraction_thresh.  (default: 0.10)
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


# Find all the fastq files in the dataset and write the full path to the file.
echo "find ${fastq_input_dir} -name \"*${fastq_file_ext}\" -type f | sed \"s/${read1_suffix}\|${read2_suffix}//g\" | rev | cut -d '/' -f1 | rev | sort -V | uniq > ${fastq_list_file}"
find ${fastq_input_dir} -name "*${fastq_file_ext}" -type f | sed "s/${read1_suffix}\|${read2_suffix}//g" | rev | cut -d '/' -f1 | rev | sort -V | uniq > ${fastq_list_file}


#step 1: trimmomatic

# Activate the trimmomatic conda environment.
conda activate trimmomatic_env

# Make the trimmomatic trim output directory.
trim_fastq_dir="${preprocessing_output_dir}/trim"
mkdir -p $trim_fastq_dir

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";
    
    fastq_read1_file="${fastq_input_dir}/${fastq_filename}${read1_suffix}";
    fastq_read2_file="${fastq_input_dir}/${fastq_filename}${read2_suffix}";

    trim_paired_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.fq";
    trim_unpaired_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.fq";
    trim_paired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.reverse.fq";
    trim_unpaired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.reverse.fq";
    
    if [ ! -s $trim_paired_fastq_file ] && [ ! -s $trim_unpaired_fastq_file ] && [  ! -s $trim_paired_rev_fastq_file ] && [  ! -s $trim_unpaired_rev_fastq_file ];
    then
        echo "trimmomatic PE ${fastq_read1_file} ${fastq_read2_file} ${trim_paired_fastq_file} ${trim_unpaired_fastq_file} ${trim_paired_rev_fastq_file} ${trim_unpaired_rev_fastq_file} LEADING:${cut_leading_bases} TRAILING:${cut_trailing_bases} SLIDINGWINDOW:${sliding_window_size}:${sliding_window_phred_score} MINLEN:${min_read_length}";
        
        trimmomatic PE ${fastq_read1_file} ${fastq_read2_file} ${trim_paired_fastq_file} ${trim_unpaired_fastq_file} ${trim_paired_rev_fastq_file} ${trim_unpaired_rev_fastq_file} LEADING:${cut_leading_bases} TRAILING:${cut_trailing_bases} SLIDINGWINDOW:${sliding_window_size}:${sliding_window_phred_score} MINLEN:${min_read_length};
    else
    
        trim_paired_fastq_file_filename=$(basename $trim_paired_fastq_file)
        trim_unpaired_fastq_file_filename=$(basename $trim_unpaired_fastq_file)
        trim_paired_rev_fastq_file_filename=$(basename $trim_paired_rev_fastq_file)
        trim_unpaired_rev_fastq_file_filename=$(basename $trim_unpaired_rev_fastq_file)
        
        echo "The following files have already been created."
        echo "${trim_paired_fastq_file_filename}"
        echo "${trim_unpaired_fastq_file_filename}"
        echo "${trim_paired_rev_fastq_file_filename}"
        echo "${trim_unpaired_rev_fastq_file_filename}"
        echo "Skipping to next set of commands!!!"
    fi
done


# Make the insert size output directory.
insert_size_dir="${output_dir}/insert_size"
mkdir -p $insert_size_dir

for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";
    
    fastq_read1_file="${fastq_input_dir}/${fastq_filename}${read1_suffix}";
    fastq_read2_file="${fastq_input_dir}/${fastq_filename}${read2_suffix}";

    trim_paired_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.fq";
    trim_unpaired_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.fq";
    trim_paired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.reverse.fq";
    trim_unpaired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.reverse.fq";
    
    insert_size_sam_file="${insert_size_dir}/${fastq_filename}_insert_size.sam"
    insert_size_bam_file="${insert_size_dir}/${fastq_filename}_insert_size.bam"
    insert_size_stats_file="${insert_size_dir}/${fastq_filename}_insert_size_stats.txt"
    
    conda activate bowtie2_env
    
    if [ ! -s $insert_size_sam_file ];
    then
        echo "bowtie2 -x ${emirge_bowtie_db} -p ${num_threads} -1 ${trim_paired_fastq_file} -2 ${trim_paired_rev_fastq_file} > ${insert_size_sam_file}"
        bowtie2 -x ${emirge_bowtie_db} -p ${num_threads} -1 ${trim_paired_fastq_file} -2 ${trim_paired_rev_fastq_file} > ${insert_size_sam_file}
    else
        insert_size_sam_file_filename=$(basename $insert_size_sam_file)
        echo "The ${insert_size_sam_file_filename} file has already been created. Skipping to next set of commands!!!"
    fi
    
    conda activate samtools_env
    
    if [ ! -s $insert_size_bam_file ];
    then
        echo "samtools view -bS < ${insert_size_sam_file} > ${insert_size_bam_file}"
        samtools view -bS < ${insert_size_sam_file} > ${insert_size_bam_file}
    else
        insert_size_bam_file_filename=$(basename $insert_size_sam_file)
        echo "The ${insert_size_bam_file_filename} file has already been created. Skipping to next set of commands!!!"
    fi
    
    if [ ! -s $insert_size_stats_file ];
    then
        echo "samtools stats ${insert_size_bam_file} > ${insert_size_stats_file}"
        samtools stats ${insert_size_bam_file} > ${insert_size_stats_file}
    else
        insert_size_stats_file_filename=$(basename $insert_size_stats_file)
        echo "The ${insert_size_stats_file_filename} file has already been created. Skipping to next set of commands!!!"
    fi
    
done




for fastq_filename in $(cat $fastq_list_file);
do
    echo "Processing ${fastq_filename}......";
    
    trim_paired_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.fq";
    trim_unpaired_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.fq";
    trim_paired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.paired.reverse.fq";
    trim_unpaired_rev_fastq_file="${trim_fastq_dir}/${fastq_filename}.unpaired.reverse.fq";
        
    insert_size_stats_file="${insert_size_dir}/${fastq_filename}_insert_size_stats.txt"
    
    # MAX_READ_LENGTH
    # length of longest read in input data.
    #SN      maximum length: 151
    max_read_length=$(grep "^SN\s\+maximum length:" ${insert_size_stats_file} | sed -e 's/SN\tmaximum length:\t\([0-9]\+\)/\1/g')
    
    # INSERT_MEAN
    # insert size distribution mean.
    #SN      insert size average:    227.2
    insert_mean=$(grep "^SN\s\+insert size average:" ${insert_size_stats_file} | sed -e 's/SN\tinsert size average:\t\([0-9]\+\)/\1/g' | cut -d '.' -f1)
    
    # INSERT_STDDEV
    # insert size distribution standard deviation.
    #SN      insert size standard deviation: 79.0
    insert_stddev=$(grep "^SN\s\+insert size standard deviation:" ${insert_size_stats_file} | sed -e 's/SN\tinsert size standard deviation:\t\([0-9]\+\)/\1/g' | cut -d '.' -f1)
        
    emirge_sample_output_dir="${emirge_output_dir}/${fastq_filename}"
    mkdir -p ${emirge_sample_output_dir}
    
    conda activate emirge_env
    
    echo "emirge.py ${emirge_sample_output_dir} -1 ${trim_paired_fastq_file} -2 ${trim_paired_rev_fastq_file} -f ${emirge_fasta_db} -b  ${emirge_bowtie_db} -l ${max_read_length} -i ${insert_mean} -s ${insert_stddev} -n ${num_iter} -a ${num_threads} -p  ${snp_fraction_thresh} -v ${variant_fraction_thresh} -j ${join_threshold} -c ${min_depth} --phred33"
    emirge.py ${emirge_sample_output_dir} -1 ${trim_paired_fastq_file} -2 ${trim_paired_rev_fastq_file} -f ${emirge_fasta_db} -b  ${emirge_bowtie_db} -l ${max_read_length} -i ${insert_mean} -s ${insert_stddev} -n ${num_iter} -a ${num_threads} -p  ${snp_fraction_thresh} -v ${variant_fraction_thresh} -j ${join_threshold} -c ${min_depth} --phred33

    emirge_iter_file="${emirge_sample_output_dir}/${fastq_filename}/iter.${num_iter}"
    emirge_renamed_file="${emirge_sample_output_dir}/${fastq_filename}/${fastq_filename}_emirge.fasta"
    
    echo "emirge_rename_fasta.py ${emirge_iter_file} > ${emirge_renamed_file}"
    emirge_rename_fasta.py ${emirge_iter_file} > ${emirge_renamed_file}
    
done


