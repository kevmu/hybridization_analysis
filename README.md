# hybridization_analysis

# Clone the hybridization_analysis github repository.
cd $HOME
git clone https://github.com/kevmu/hybridization_analysis.git

cd hybridization_analysis

# Install Conda environments.

conda create --name cutadapt_env
conda activate cutadapt_env
conda install -c bioconda cutadapt

conda create --name biopython_env
conda activate biopython_env
conda install -c bioconda biopython

conda create --name cd_hit_env
conda activate cd_hit_env
conda install -c bioconda cd-hit

conda create --name cutadapt_env
conda activate cutadapt_env
conda install -c bioconda cutadapt

conda create --name prinseq_env
conda activate prinseq_env
conda install -c bioconda prinseq

conda create --name emirge_env
conda activate emirge_env
conda install -c bioconda emirge

conda create --name vsearch_env
conda activate vsearch_env
conda install -c bioconda vsearch

conda create --name qiime_env
conda activate qiime_env
conda install -c conda-forge -c bioconda qiime

## Download source code.

mkdir -p software
cd software

## EMIRGE
# Clone the EMIRGE source code.
git clone https://github.com/csmiller/EMIRGE.git

# Download the prinseq source code.
wget https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download -O prinseq-lite-0.20.4.tar.gz

## PRINSEQ
# Uncompress the prinseq-lite-0.20.4.tar.gz tar gzipped file.
tar xvzf prinseq-lite-0.20.4.tar.gz

# Download the usearch source code.
wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz -O usearch.gz

# Add the path to usearch to your ~/.bashrc file.
export PATH="$HOME/hybridization_analysis/software/prinseq-lite-0.20.4:$PATH"

## usearch

# Uncompress the usearch binary code.
gunzip -c usearch.gz > usearch

# Make usearch permissions all read and execute.
chmod a+rx usearch

# Add the path to usearch to your ~/.bashrc file.
export PATH="$HOME/hybridization_analysis/software:$PATH"

# Get the source of your .bashrc or .bash_profile
source ~/.bashrc

# OR
source ~/.bash_profile

# Install EMIRGE 16S SILVA database.
sh create_emirge_databases.sh &> run_create_emirge_databases.log.txt

# Get the filenames of each sample so that we can generate a fasta file list of sample names without the R1 and R2 fastq read suffix.

# Example using the find command.
find /archive/dumonceauxt/230801_M01666_0203_000000000-KTT65_IMPACTTlilacBBSPTWMQhybsEPL/Fastq -type f -name "*.fastq.gz" | sed s/_L001_R1_001.fastq.gz|_L001_R2_001.fastq.gz//g | rev | cut -d '/' -f1 | rev | sort -V | uniq | grep 16S > /home/AGR.GC.CA/muirheadk/hybridization_analysis/fastq_files_list.txt

