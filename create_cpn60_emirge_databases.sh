#!/bin/bash
#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=60G
#SBATCH --output=run_create_emirge_databases.%A.out
#SBATCH --error=run_create_emirge_databases.%A.err

source ~/.bashrc

# The cpn60 input fasta file.
cpn60_db_fasta="/home/AGR.GC.CA/muirheadk/hybridization_analysis/cpn60_all_nut_seq.txt"

db_dir="/home/AGR.GC.CA/muirheadk/hybridization_analysis/emirge_dbs"

emirge_scripts_dir="/home/AGR.GC.CA/muirheadk/hybridization_analysis/software/EMIRGE"

# The number of threads to use in cd-hit.
num_threads=8

cpn60_db_dir="${db_dir}/cpn60_db"
mkdir -p $cpn60_db_dir

# Activate the cd-hit conda environment.
conda activate cd_hit_env

clustered_cpn60_db_fasta="${cpn60_db_dir}/cpn60_all_nut_seq.clustered.emirge.ref.fasta"

if [ ! -s $clustered_cpn60_db_fasta ];
then
	echo -e "cd-hit -i ${cpn60_db_fasta} -c 1.0 -d 3000 -M 5000 -T ${num_threads} -o ${clustered_cpn60_db_fasta}\n"

	cd-hit -i ${cpn60_db_fasta} -c 1.0 -d 3000 -M 5000 -T ${num_threads} -o ${clustered_cpn60_db_fasta}
else
	clustered_cpn60_db_fasta_filename=$(basename $clustered_cpn60_db_fasta)
	echo "The ${clustered_cpn60_db_fasta_filename} file has already been created. Skipping to next set of commands!!!"
fi

# Activate the emirge conda environment.
conda activate emirge_env


emirge_cpn60_db_prefix="${cpn60_db_dir}/cpn60_all_nut_seq.clustered.emirge.ref"
##find /bulk/sycuro_bulk/lsycuro_labshare/kevin/testing/hybridization_analysis/dbs/emirge_db/cpn60_db -type f -name "*.bt2"

# I NEED TO ADD THE IF STATEMENT SFOR THE BOWTIE DATABASE EXTENSION FILES
#if [ ! -s $emirge_cpn60_db_fasta ];
#then

echo "bowtie-build ${clustered_cpn60_db_fasta} ${emirge_cpn60_db_prefix}"
bowtie-build ${clustered_cpn60_db_fasta} ${emirge_cpn60_db_prefix}

#else
#	emirge_cpn60_db_fasta_filename=$(basename $emirge_cpn60_db_fasta)
#	echo "The ${_filename} file has already been created. Skipping to next set of commands!!!"
#fi


echo "The emirge databases have been created successfully!!! ${emirge_cpn60_db_fasta} and ${emirge_cpn60_db_prefix}."


