#Clustering of similar protein sequences
conda activate mmseqs2-14.7e284
mmseqs createdb Prot_sequences.fasta MMSEQ_DB
mmseqs linclust MMSEQ_DB DB_clu tmp --min-seq-id 0.99 --cov-mode 1 -c 0.8
mmseqs createsubdb DB_clu MMSEQ_DB DB_clu_rep
mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta

#Alignment of protein sequences
conda activate mafft-7.487
mafft --auto DB_clu_rep.fasta > Prot.aln
conda deactivate
##Manual removal of false postive sequences

#Addittion of reference sequences (siderophore NRPS or NIS protein sequences)
cat DB_clu_rep_filt.fasta references.fasta > Prot_and_references.fasta

#Alignment of protein sequences and references
conda activate mafft-7.487
mafft --auto Prot_and_references.fasta > Prot_and_references.aln
conda deactivate

#Construction of a phylogenetic tree based on the alignment of protein sequences
conda activate iqtree-2.2.0.3
iqtree2 -s Prot_and_references.aln --seqtype AA -T 20 -m TEST -bb 1000 -alrt 1000
conda deactivate
##Manual removal of false postive sequences

#Alignment of filtered protein sequences
conda activate mafft-7.487
mafft --auto Prot_final.fasta > Prot_final.aln
conda deactivate

#Contruction of HMM, construct binary compressed datafiles for hmmscan, run of hmmscan on control sequences and metagenomes
conda activate hmmer-3.3.2
hmmbuild X.hmm Prot_final.aln
hmmpress X.hmm
hmmscan --domtblout Control.tbl X.hmm Prot_and_references.fasta
for file in /METAGENOMES/AOP*.faa;
do

   samplename=$(basename $file .faa)
   hmmscan --domtblout ${samplename}.tbl -E 0.00001 X.hmm /METAGENOMES/${samplename}.faa

done
conda deactivate
##Set of bitscore below last positive hit