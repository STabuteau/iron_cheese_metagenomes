This folder contains 9 files:
# bash scripts
HMMbuilt.sh is a script used to design and test HMM
HMMer_MAGs_Genomes.sh is a script used to search HMMs on MAGs and/or genomes using HMMer
HMMer_Metag.sh is a script used to search HMMs on metagenomes using HMMer and normalise data in RPKM
fungal_MAGs_analysis.sh is a script used to identify fungal MAGs and to classify them taxonomically

# metadata file
HMMs_metadata.tsv is a metadata file containing information about the HMMs `iron_cheese_metagenomes/HMMs/`

# R scripts
append_cov.R is a R script used in HMMer_Metag.sh
merge_cov_qlen.R is a R script used in HMMer_Metag.sh
merge_files.R is a R script used in HMMer_Metag.sh HMMer_MAGs_Genomes.sh fungal_MAGs_analysis.sh

# perl script
catfasta2phyml.pl is a script used in fungal_MAGs_analysis.sh to concatenate gene marker alignment, this script belong to Johan Nylander (https://github.com/nylander/catfasta2phyml)
