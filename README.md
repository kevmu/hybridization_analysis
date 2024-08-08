# hybridization_analysis
\
conda create --name emirge_env\
conda activate emirge_env\
conda install -c bioconda emirge\
\
conda create --name trimmomatic_env\
conda activate trimmomatic_env\
conda install -c bioconda trimmomatic\
\
conda create --name biopython_env\
conda activate biopython_env\
conda install -c bioconda biopython_env\
\
conda create --name cd_hit_env\
conda activate cd_hit_env\
conda install -c bioconda cd-hit\

conda create --name cutadapt_env
conda activate cutadapt_env
conda install -c bioconda cutadapt

conda create --name prinseq_env
conda activate prinseq_env
conda install -c bioconda prinseq

tar xvzf prinseq-lite-0.20.4.tar.gz
wget https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download -O prinseq-lite-0.20.4.tar.gz

wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
cp usearch11.0.667_i86linux32 usearch

# Make usearch permissions all read and execute.
chmod a+rx usearch

# Add the path to usearch to your .bashrc or .bash_profile
export PATH=/work/sycuro_lab/kevin/software:$PATH


# Get the source of your .bashrc or .bash_profile
source ~/.bashrc

# OR
source ~/.bash_profile

