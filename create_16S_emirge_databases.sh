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

db_dir="$HOME/hybridization_analysis/databases"

emirge_db_dir="${db_dir}/emirge_dbs"

emirge_scripts_dir="$HOME/hybridization_analysis/software/EMIRGE"

# The number of threads to use in cd-hit.
num_threads=8

emirge_silva_db_dir="${emirge_db_dir}/silva_db"
mkdir -p $emirge_silva_db_dir

compressed_silva_db="${emirge_silva_db_dir}/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz"
silva_db_fasta="${emirge_silva_db_dir}/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta"

if [ ! -s $compressed_silva_db ];
then

	echo -e "wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz -O ${compressed_silva_db}\n"
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta.gz -O ${compressed_silva_db}
else
	compressed_silva_db_filename=$(basename $compressed_silva_db)
	echo "The ${compressed_silva_db_filename} file has already been created. Skipping to next set of commands!!!"
fi

if [ ! -s $silva_db_fasta ];
then

	echo -e "gunzip -c ${compressed_silva_db} > ${silva_db_fasta}\n"
	gunzip -c ${compressed_silva_db} > ${silva_db_fasta}
else
	silva_db_fasta_filename=$(basename $silva_db_fasta)
	echo "The ${silva_db_fasta_filename} file has already been created. Skipping to next set of commands!!!"
fi

# Activate the emirge conda environment.
conda activate emirge_env

silva_fixed_db_fasta="${emirge_silva_db_dir}/SILVA_138_SSURef_NR99_tax_silva_trunc.fixed.fasta"

if [ ! -s $silva_fixed_db_fasta ];
then

	# Fix the nonstandard characters in the silva database.
	echo "python ${emirge_scripts_dir}/utils/fix_nonstandard_chars.py < ${silva_db_fasta} > ${silva_fixed_db_fasta}"
	python ${emirge_scripts_dir}/utils/fix_nonstandard_chars.py < ${silva_db_fasta} > ${silva_fixed_db_fasta}
else
	silva_fixed_db_fasta_filename=$(basename $silva_fixed_db_fasta)
	echo "The ${silva_fixed_db_fasta_filename} file has already been created. Skipping to next set of commands!!!"
fi

# Activate the cd-hit conda environment.
conda activate cd_hit_env

clustered_silva_db_fasta="${emirge_silva_db_dir}/SILVA_138_SSURef_NR97_tax_silva_trunc.fixed.clustered.fasta"

if [ ! -s $clustered_silva_db_fasta ];
then
	echo -e "cd-hit -i ${silva_fixed_db_fasta} -c 0.97 -d 3000 -M 5000 -T ${num_threads} -o ${clustered_silva_db_fasta}\n"

	cd-hit -i ${silva_fixed_db_fasta} -c 0.97 -d 3000 -M 5000 -T ${num_threads} -o ${clustered_silva_db_fasta}
else
	clustered_silva_db_fasta_filename=$(basename $clustered_silva_db_fasta)
	echo "The ${clustered_silva_db_fasta_filename} file has already been created. Skipping to next set of commands!!!"
fi

# Activate the biopython conda environment.
conda activate biopython_env

emirge_silva_db_fasta="${emirge_silva_db_dir}/SILVA_138_SSURef_NR97_tax_silva_trunc.fixed.clustered.emirge.ref.fasta"

if [ ! -s $emirge_silva_db_fasta ];
then

	echo "python filter_silva_emirge_db_seqs.py --fasta_infile ${clustered_silva_db_fasta} --output_dir ${emirge_silva_db_dir}"
	python filter_silva_emirge_db_seqs.py --fasta_infile ${clustered_silva_db_fasta} --output_dir ${emirge_silva_db_dir}

else
	emirge_silva_db_fasta_filename=$(basename $emirge_silva_db_fasta)
	echo "The ${emirge_silva_db_fasta_filename} file has already been created. Skipping to next set of commands!!!"
fi


# Activate the emirge conda environment.
conda activate emirge_env

emirge_silva_db_prefix="${emirge_silva_db_dir}/SILVA_138_SSURef_NR97_tax_silva_trunc.fixed.clustered.emirge.ref"
##find /bulk/sycuro_bulk/lsycuro_labshare/kevin/testing/hybridization_analysis/dbs/emirge_db/silva_db -type f -name "*.bt2"

# I NEED TO ADD THE IF STATEMENT SFOR THE BOWTIE DATABASE EXTENSION FILES
#if [ ! -s $emirge_silva_db_fasta ];
#then

echo "bowtie-build ${emirge_silva_db_fasta} ${emirge_silva_db_prefix}"
bowtie-build ${emirge_silva_db_fasta} ${emirge_silva_db_prefix}

#else
#	emirge_silva_db_fasta_filename=$(basename $emirge_silva_db_fasta)
#	echo "The ${_filename} file has already been created. Skipping to next set of commands!!!"
#fi

bowtie2_db_dir="${db_dir}/bowtie2_dbs"
mkdir -p $bowtie2_db_dir

bowtie2_silva_db_dir="${bowtie2_db_dir}/silva_db"
mkdir -p $bowtie2_silva_db_dir

bowtie2_silva_db_fasta="${bowtie2_silva_db_dir}/SILVA_138_SSURef_NR97_tax_silva_trunc.fixed.clustered.emirge.ref.fasta"
bowtie2_silva_db_prefix="${bowtie2_silva_db_dir}/SILVA_138_SSURef_NR97_tax_silva_trunc.fixed.clustered.emirge.ref"

echo "cp ${emirge_silva_db_fasta} ${bowtie2_silva_db_fasta}"
cp ${emirge_silva_db_fasta} ${bowtie2_silva_db_fasta}

# Activate the bowtie2 conda environment.
conda activate bowtie2_env

# I NEED TO ADD THE IF STATEMENT S^?FOR THE BOWTIE DATABASE EXTENSION FILES
#if [ ! -s $emirge_silva_db_fasta ];
#then

echo "bowtie2-build ${bowtie2_silva_db_fasta} ${bowtie2_silva_db_prefix}" 
bowtie2-build ${bowtie2_silva_db_fasta} ${bowtie2_silva_db_prefix}

#else
#       emirge_silva_db_fasta_filename=$(basename $emirge_silva_db_fasta)
#       echo "The ${_filename} file has already been created. Skipping to next set of commands!!!"
#fi

vsearch_gtdb_dir="${db_dir}/vsearch_gtdb"
mkdir -p $vsearch_gtdb_dir

compressed_vsearch_gtdb="${vsearch_gtdb_dir}/bac120_ssu_reps_r220.fna.gz"
vsearch_gtdb_fasta="${vsearch_gtdb_dir}/bac120_ssu_reps_r220.fna"

echo "wget \"https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/bac120_ssu_reps_r220.fna.gz\" -O ${compressed_vsearch_gtdb}"
wget "https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_reps/bac120_ssu_reps_r220.fna.gz" -O ${compressed_vsearch_gtdb}

echo "gunzip -c ${compressed_vsearch_gtdb} > ${vsearch_gtdb_fasta}"
gunzip -c ${compressed_vsearch_gtdb} > ${vsearch_gtdb_fasta}

# Activate the cd-hit conda environment.
conda activate cd_hit_env

clustered_vsearch_gtdb_fasta="${vsearch_gtdb_dir}/gtdb_bac120_ssu_reps_r220_vsearch.fasta"

echo "cd-hit -i ${vsearch_gtdb_fasta} -c 1.0 -d 3000 -M 5000 -T ${num_threads} -o ${clustered_vsearch_gtdb_fasta}"
cd-hit -i ${vsearch_gtdb_fasta} -c 1.0 -d 3000 -M 5000 -T ${num_threads} -o ${clustered_vsearch_gtdb_fasta}

echo "The emirge databases have been created successfully!!! ${emirge_silva_db_fasta} and ${emirge_silva_db_prefix}."
echo "The bowtie2 databases have been created successfully!!! ${bowtie2_silva_db_fasta} and ${bowtie2_silva_db_prefix}."
echo "The bowtie2 databases have been created successfully!!! ${clustered_vsearch_gtdb_fasta}"




