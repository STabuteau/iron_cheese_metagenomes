This folder contains 9 files:
# bash scripts
HMMbuilt.sh is a script used to design and test HMM.  
HMMer_MAGs_Genomes.sh is a script used to search HMMs on MAGs and/or genomes using HMMer.  
HMMer_Metag.sh is a script used to search HMMs on metagenomes using HMMer and normalise data in RPKM.  
fungal_MAGs_analysis.sh is a script used to identify fungal MAGs and to classify them taxonomically.  

# metadata file
HMMs_metadata.tsv is a metadata file containing information about the HMMs `iron_cheese_metagenomes/HMMs/`.  

# R scripts (v4.4.1)
append_cov.R is a R script used in HMMer_Metag.sh to add gene abundance to HMMer output file.  
merge_cov_qlen.R is a R script used in HMMer_Metag.sh to add mean querry lenght to HMMer output file.  
merge_files.R is a R script used in HMMer_Metag.sh ; HMMer_MAGs_Genomes.sh and fungal_MAGs_analysis.sh to add HMM metadata to HMMer output file.  

# perl script
catfasta2phyml.pl is a script used in fungal_MAGs_analysis.sh to concatenate gene marker alignment, this script belongs to Johan Nylander (https://github.com/nylander/catfasta2phyml).  
