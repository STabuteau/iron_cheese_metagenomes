#Identification of fungal MAGs
conda activate busco-5.3.2
for file in /Potential_fungal_fasta/AOP*.fa;
do
   samplename=$(basename $file .fa)
   busco -i /Potential_fungal_fasta/${samplename}.fa --auto-lineage-euk -o ${samplename} -m genome
done

#Annotation of fungal MAGs
conda activate metaeuk-6.a5d39d9
for file in /Fungal_fasta/AOP*.fa;
do
   samplename=$(basename $file .fa)
   metaeuk easy-predict --metaeuk-eval 0.0001 --metaeuk-tcov 0.6 --min-length 20 /Fungal_fasta/${samplename}.fa /UniRef100_db/UniRef100 METAEUK/${samplename} ${samplename}_tempFolder
   rm -fr ${samplename}_tempFolder
done

#Construct binary compressed datafiles for hmmscan, run of hmmscan and filter of results by predefined bitscore cutoffs for each HMM and length of the CDS recognised by the HMMs
conda activate hmmer-3.3.2
hmmpress All_HMMv4_fungal_ver.hmm
for file in /METAEUK/AOP*.fas;
do
   samplename=$(basename $file .fas)
   hmmscan --domtblout /HMMER/${samplename}.tbl -E 0.00001 All_HMMv4_fungal_ver.hmm /METAEUK/${samplename}.fas
done
conda deactivate
for file in /HMMER/AOP*.tbl;
do
   samplename=$(basename $file .tbl)
   cat ${samplename}.tbl | grep -v '^#' | \
   sed 's/\s\s*/\t/g' | \
   sort -k4,4 | \
   awk 'BEGIN{FS="\t"; OFS="\t"}{if(($19-$18+1)/$6>=0.7){print $0,($19-$18+1),($19-$18+1)/$6}}' | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,$3,$4,$6,$8,$13,$16,$17,$18,$19,$NF}' | \
   sed 's/,/\./g' > ${samplename}_parse_file.tmp
   Rscript merge_files.R ${samplename}_parse_file.tmp HMMs_metadata.tsv ${samplename}_parse_merge.tsv
   awk -F '\t' '{if($5>=$NF) {print $0}}' ${samplename}_parse_merge.tsv | awk '$5 > max[$3] { max[$3] = $5; m[$3] = $0 } END { for (i in m) { print m[i] } }' > ${samplename}_parse_merge_bitscore.tsv
done
rm -f *_parse_file.tmp
rm -f *_parse_merge.tsv

#Identification of marker genes in MAGs and supplemantary references
conda activate ufcg-1.0.6
for file in /Fungal_fasta_and_references/*.fa;
do
   samplename=$(basename $file .fa)
   ufcg profile -i /Fungal_fasta_and_references/${samplename}.fa -o /UFCG/${samplename}
done
conda deactivate

#Alignment of gene marker of MAGs, supplementary references and UFCG references
conda activate ufcg-1.0.6
ufcg align -i UFCG/ -o UFCG/align/
conda deactivate
conda activate mafft-7.487
for file in /UFCG/align/*.fasta;
do
   samplename=$(basename $file .zZ.fasta)
   mafft --auto ${samplename}.zZ.fasta > ${samplename}.aln
done
conda deactivate

#Trimming of each alignment
conda activate trimal-1.4.1
for file in /UFCG/align/*_pro.aln;
do
   samplename=$(basename $file _pro.aln)
   trimal -in ${samplename}_pro.aln -out ${samplename}_pro_trim.aln -gt 0.1
done
conda deactivate

#Concatenation of gene marker alignment
##the perl script was obtained from https://github.com/nylander/catfasta2phyml
perl catfasta2phyml.pl --concatenate -f align/*_pro_trim.aln > align/all_MAGs_pro_trim.aln


#Construction of a phylogenetic tree based on the alignment of marker genes
conda activate iqtree-2.2.0.3 
iqtree2 -s /UFCG/align/all_MAGS_pro_trim.aln --seqtype AA -T 36 -m TEST -bb 1000 -alrt 1000
conda deactivate

#Calculation of the relative abundance of fungal MAGs
conda activate coverm-0.6.1
coverm genome --coupled /Raw_Reads/AOP10_L1_surf_R1.fastq.gz /Raw_Reads/AOP10_L1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP10_L1_surf_bin.4.fa -o AOP10_L1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP11_C1_surf_R1.fastq.gz /Raw_Reads/AOP11_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP11_C1_surf_bin.6.fa -o AOP11_C1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP11_C1_surf_R1.fastq.gz /Raw_Reads/AOP11_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP11_C1_surf_bin.9.fa -o AOP11_C1_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP12_E2_surf_R1.fastq.gz /Raw_Reads/AOP12_E2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP12_E2_surf_bin.4.fa -o AOP12_E2_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP12_E2_surf_R1.fastq.gz /Raw_Reads/AOP12_E2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP12_E2_surf_bin.5.fa -o AOP12_E2_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP13_A1_surf_R1.fastq.gz /Raw_Reads/AOP13_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP13_A1_surf_bin.4.fa -o AOP13_A1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP13_A1_surf_R1.fastq.gz /Raw_Reads/AOP13_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP13_A1_surf_bin.5.fa -o AOP13_A1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP13_A3_surf_R1.fastq.gz /Raw_Reads/AOP13_A3_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP13_A3_surf_bin.2.fa -o AOP13_A3_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP13_A3_surf_R1.fastq.gz /Raw_Reads/AOP13_A3_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP13_A3_surf_bin.6.fa -o AOP13_A3_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP14_B2_surf_R1.fastq.gz /Raw_Reads/AOP14_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP14_B2_surf_bin.1.fa -o AOP14_B2_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP14_B2_surf_R1.fastq.gz /Raw_Reads/AOP14_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP14_B2_surf_bin.2.fa -o AOP14_B2_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP14_E2_surf_R1.fastq.gz /Raw_Reads/AOP14_E2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP14_E2_surf_bin.12.fa -o AOP14_E2_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP15_B1_surf_R1.fastq.gz /Raw_Reads/AOP15_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP15_B1_surf_bin.6.fa -o AOP15_B1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP15_C2_surf_R1.fastq.gz /Raw_Reads/AOP15_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP15_C2_surf_bin.8.fa -o AOP15_C2_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP15_D2_surf_R1.fastq.gz /Raw_Reads/AOP15_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP15_D2_surf_bin.2.fa -o AOP15_D2_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP15_D2_surf_R1.fastq.gz /Raw_Reads/AOP15_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP15_D2_surf_bin.4.fa -o AOP15_D2_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP16_AQ2_surf_R1.fastq.gz /Raw_Reads/AOP16_AQ2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP16_AQ2_surf_bin.1.fa -o AOP16_AQ2_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP17_E1_surf_R1.fastq.gz /Raw_Reads/AOP17_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP17_E1_surf_bin.10.fa -o AOP17_E1_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP17_E1_surf_R1.fastq.gz /Raw_Reads/AOP17_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP17_E1_surf_bin.6.fa -o AOP17_E1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP17_E1_surf_R1.fastq.gz /Raw_Reads/AOP17_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP17_E1_surf_bin.9.fa -o AOP17_E1_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP18_A1_surf_R1.fastq.gz /Raw_Reads/AOP18_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_A1_surf_bin.11.fa -o AOP18_A1_surf_bin.11
coverm genome --coupled /Raw_Reads/AOP18_A1_surf_R1.fastq.gz /Raw_Reads/AOP18_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_A1_surf_bin.9.fa -o AOP18_A1_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP18_A2_surf_R1.fastq.gz /Raw_Reads/AOP18_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_A2_surf_bin.22.fa -o AOP18_A2_surf_bin.22
coverm genome --coupled /Raw_Reads/AOP18_A2_surf_R1.fastq.gz /Raw_Reads/AOP18_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_A2_surf_bin.8.fa -o AOP18_A2_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP18_B2_surf_R1.fastq.gz /Raw_Reads/AOP18_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_B2_surf_bin.10.fa -o AOP18_B2_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP18_B2_surf_R1.fastq.gz /Raw_Reads/AOP18_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_B2_surf_bin.13.fa -o AOP18_B2_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP18_B2_surf_R1.fastq.gz /Raw_Reads/AOP18_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_B2_surf_bin.14.fa -o AOP18_B2_surf_bin.14
coverm genome --coupled /Raw_Reads/AOP18_B2_surf_R1.fastq.gz /Raw_Reads/AOP18_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_B2_surf_bin.8.fa -o AOP18_B2_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP18_D1_surf_R1.fastq.gz /Raw_Reads/AOP18_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_D1_surf_bin.16.fa -o AOP18_D1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP18_D1_surf_R1.fastq.gz /Raw_Reads/AOP18_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_D1_surf_bin.4.fa -o AOP18_D1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP18_D1_surf_R1.fastq.gz /Raw_Reads/AOP18_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP18_D1_surf_bin.9.fa -o AOP18_D1_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP19_D1_surf_R1.fastq.gz /Raw_Reads/AOP19_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP19_D1_surf_bin.8.fa -o AOP19_D1_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP19_E1_surf_R1.fastq.gz /Raw_Reads/AOP19_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP19_E1_surf_bin.16.fa -o AOP19_E1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP19_H1_surf_R1.fastq.gz /Raw_Reads/AOP19_H1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP19_H1_surf_bin.1.fa -o AOP19_H1_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP19_H1_surf_R1.fastq.gz /Raw_Reads/AOP19_H1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP19_H1_surf_bin.4.fa -o AOP19_H1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP1_A1_surf_R1.fastq.gz /Raw_Reads/AOP1_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_A1_surf_bin.20.fa -o AOP1_A1_surf_bin.20
coverm genome --coupled /Raw_Reads/AOP1_A1_surf_R1.fastq.gz /Raw_Reads/AOP1_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_A1_surf_bin.5.fa -o AOP1_A1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP1_C1_surf_R1.fastq.gz /Raw_Reads/AOP1_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_C1_surf_bin.16.fa -o AOP1_C1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP1_E1_surf_R1.fastq.gz /Raw_Reads/AOP1_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_E1_surf_bin.32.fa -o AOP1_E1_surf_bin.32
coverm genome --coupled /Raw_Reads/AOP1_G1_surf_R1.fastq.gz /Raw_Reads/AOP1_G1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_G1_surf_bin.13.fa -o AOP1_G1_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP1_H1_surf_R1.fastq.gz /Raw_Reads/AOP1_H1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_H1_surf_bin.5.fa -o AOP1_H1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP1_I1_surf_R1.fastq.gz /Raw_Reads/AOP1_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP1_I1_surf_bin.18.fa -o AOP1_I1_surf_bin.18
coverm genome --coupled /Raw_Reads/AOP20_C1_surf_R1.fastq.gz /Raw_Reads/AOP20_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP20_C1_surf_bin.6.fa -o AOP20_C1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP21_A1_surf_R1.fastq.gz /Raw_Reads/AOP21_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP21_A1_surf_bin.1.fa -o AOP21_A1_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP21_A1_surf_R1.fastq.gz /Raw_Reads/AOP21_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP21_A1_surf_bin.6.fa -o AOP21_A1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP21_I1_surf_R1.fastq.gz /Raw_Reads/AOP21_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP21_I1_surf_bin.3.fa -o AOP21_I1_surf_bin.3
coverm genome --coupled /Raw_Reads/AOP22_C1_surf_R1.fastq.gz /Raw_Reads/AOP22_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP22_C1_surf_bin.7.fa -o AOP22_C1_surf_bin.7
coverm genome --coupled /Raw_Reads/AOP22_E1_surf_R1.fastq.gz /Raw_Reads/AOP22_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP22_E1_surf_bin.12.fa -o AOP22_E1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP22_E1_surf_R1.fastq.gz /Raw_Reads/AOP22_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP22_E1_surf_bin.17.fa -o AOP22_E1_surf_bin.17
coverm genome --coupled /Raw_Reads/AOP23_C1_surf_R1.fastq.gz /Raw_Reads/AOP23_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP23_C1_surf_bin.5.fa -o AOP23_C1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP23_C1_surf_R1.fastq.gz /Raw_Reads/AOP23_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP23_C1_surf_bin.6.fa -o AOP23_C1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP23_I1_surf_R1.fastq.gz /Raw_Reads/AOP23_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP23_I1_surf_bin.10.fa -o AOP23_I1_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP23_I1_surf_R1.fastq.gz /Raw_Reads/AOP23_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP23_I1_surf_bin.14.fa -o AOP23_I1_surf_bin.14
coverm genome --coupled /Raw_Reads/AOP23_I1_surf_R1.fastq.gz /Raw_Reads/AOP23_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP23_I1_surf_bin.8.fa -o AOP23_I1_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP24_D1_surf_R1.fastq.gz /Raw_Reads/AOP24_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP24_D1_surf_bin.8.fa -o AOP24_D1_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP25_A1_surf_R1.fastq.gz /Raw_Reads/AOP25_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_A1_surf_bin.23.fa -o AOP25_A1_surf_bin.23
coverm genome --coupled /Raw_Reads/AOP25_A1_surf_R1.fastq.gz /Raw_Reads/AOP25_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_A1_surf_bin.2.fa -o AOP25_A1_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP25_A1_surf_R1.fastq.gz /Raw_Reads/AOP25_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_A1_surf_bin.6.fa -o AOP25_A1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP25_B1_coeur_R1.fastq.gz /Raw_Reads/AOP25_B1_coeur_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B1_coeur_bin.20.fa -o AOP25_B1_coeur_bin.20
coverm genome --coupled /Raw_Reads/AOP25_B1_surf_R1.fastq.gz /Raw_Reads/AOP25_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B1_surf_bin.16.fa -o AOP25_B1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP25_B2_surf_R1.fastq.gz /Raw_Reads/AOP25_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B2_surf_bin.12.fa -o AOP25_B2_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP25_B2_surf_R1.fastq.gz /Raw_Reads/AOP25_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B2_surf_bin.13.fa -o AOP25_B2_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP25_B2_surf_R1.fastq.gz /Raw_Reads/AOP25_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B2_surf_bin.16.fa -o AOP25_B2_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP25_B2_surf_R1.fastq.gz /Raw_Reads/AOP25_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B2_surf_bin.1.fa -o AOP25_B2_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP25_B2_surf_R1.fastq.gz /Raw_Reads/AOP25_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_B2_surf_bin.5.fa -o AOP25_B2_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP25_C1_surf_R1.fastq.gz /Raw_Reads/AOP25_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_C1_surf_bin.12.fa -o AOP25_C1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP25_C1_surf_R1.fastq.gz /Raw_Reads/AOP25_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_C1_surf_bin.3.fa -o AOP25_C1_surf_bin.3
coverm genome --coupled /Raw_Reads/AOP25_D2_surf_R1.fastq.gz /Raw_Reads/AOP25_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_D2_surf_bin.14.fa -o AOP25_D2_surf_bin.14
coverm genome --coupled /Raw_Reads/AOP25_D2_surf_R1.fastq.gz /Raw_Reads/AOP25_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_D2_surf_bin.18.fa -o AOP25_D2_surf_bin.18
coverm genome --coupled /Raw_Reads/AOP25_D2_surf_R1.fastq.gz /Raw_Reads/AOP25_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_D2_surf_bin.19.fa -o AOP25_D2_surf_bin.19
coverm genome --coupled /Raw_Reads/AOP25_D2_surf_R1.fastq.gz /Raw_Reads/AOP25_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_D2_surf_bin.1.fa -o AOP25_D2_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP25_D2_surf_R1.fastq.gz /Raw_Reads/AOP25_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_D2_surf_bin.8.fa -o AOP25_D2_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP25_D2_surf_R1.fastq.gz /Raw_Reads/AOP25_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_D2_surf_bin.9.fa -o AOP25_D2_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP25_F1_surf_R1.fastq.gz /Raw_Reads/AOP25_F1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_F1_surf_bin.6.fa -o AOP25_F1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP25_F2_surf_R1.fastq.gz /Raw_Reads/AOP25_F2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP25_F2_surf_bin.18.fa -o AOP25_F2_surf_bin.18
coverm genome --coupled /Raw_Reads/AOP26_D1_surf_R1.fastq.gz /Raw_Reads/AOP26_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP26_D1_surf_bin.10.fa -o AOP26_D1_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP26_D1_surf_R1.fastq.gz /Raw_Reads/AOP26_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP26_D1_surf_bin.5.fa -o AOP26_D1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP27_A1_surf_R1.fastq.gz /Raw_Reads/AOP27_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP27_A1_surf_bin.15.fa -o AOP27_A1_surf_bin.15
coverm genome --coupled /Raw_Reads/AOP28_A2_surf_R1.fastq.gz /Raw_Reads/AOP28_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP28_A2_surf_bin.19.fa -o AOP28_A2_surf_bin.19
coverm genome --coupled /Raw_Reads/AOP28_C2_surf_R1.fastq.gz /Raw_Reads/AOP28_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP28_C2_surf_bin.1.fa -o AOP28_C2_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP28_C2_surf_R1.fastq.gz /Raw_Reads/AOP28_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP28_C2_surf_bin.7.fa -o AOP28_C2_surf_bin.7
coverm genome --coupled /Raw_Reads/AOP28_E1_surf_R1.fastq.gz /Raw_Reads/AOP28_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP28_E1_surf_bin.6.fa -o AOP28_E1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP29_E1_surf_R1.fastq.gz /Raw_Reads/AOP29_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP29_E1_surf_bin.5.fa -o AOP29_E1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP2_E1_surf_R1.fastq.gz /Raw_Reads/AOP2_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP2_E1_surf_bin.1.fa -o AOP2_E1_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP30_B1_surf_R1.fastq.gz /Raw_Reads/AOP30_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP30_B1_surf_bin.12.fa -o AOP30_B1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP30_B2_surf_R1.fastq.gz /Raw_Reads/AOP30_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP30_B2_surf_bin.21.fa -o AOP30_B2_surf_bin.21
coverm genome --coupled /Raw_Reads/AOP30_C2_surf_R1.fastq.gz /Raw_Reads/AOP30_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP30_C2_surf_bin.20.fa -o AOP30_C2_surf_bin.20
coverm genome --coupled /Raw_Reads/AOP30_D2_surf_R1.fastq.gz /Raw_Reads/AOP30_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP30_D2_surf_bin.13.fa -o AOP30_D2_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP30_D2_surf_R1.fastq.gz /Raw_Reads/AOP30_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP30_D2_surf_bin.17.fa -o AOP30_D2_surf_bin.17
coverm genome --coupled /Raw_Reads/AOP30_D2_surf_R1.fastq.gz /Raw_Reads/AOP30_D2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP30_D2_surf_bin.19.fa -o AOP30_D2_surf_bin.19
coverm genome --coupled /Raw_Reads/AOP31_A2_surf_R1.fastq.gz /Raw_Reads/AOP31_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP31_A2_surf_bin.9.fa -o AOP31_A2_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP31_B2_surf_R1.fastq.gz /Raw_Reads/AOP31_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP31_B2_surf_bin.14.fa -o AOP31_B2_surf_bin.14
coverm genome --coupled /Raw_Reads/AOP33_2_surf_R1.fastq.gz /Raw_Reads/AOP33_2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP33_2_surf_bin.18.fa -o AOP33_2_surf_bin.18
coverm genome --coupled /Raw_Reads/AOP33_8_surf_R1.fastq.gz /Raw_Reads/AOP33_8_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP33_8_surf_bin.3.fa -o AOP33_8_surf_bin.3
coverm genome --coupled /Raw_Reads/AOP33_8_surf_R1.fastq.gz /Raw_Reads/AOP33_8_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP33_8_surf_bin.7.fa -o AOP33_8_surf_bin.7
coverm genome --coupled /Raw_Reads/AOP37_A2_surf_R1.fastq.gz /Raw_Reads/AOP37_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP37_A2_surf_bin.4.fa -o AOP37_A2_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP37_A2_surf_R1.fastq.gz /Raw_Reads/AOP37_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP37_A2_surf_bin.6.fa -o AOP37_A2_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP38_C1_surf_R1.fastq.gz /Raw_Reads/AOP38_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP38_C1_surf_bin.14.fa -o AOP38_C1_surf_bin.14
coverm genome --coupled /Raw_Reads/AOP39_D1_surf_R1.fastq.gz /Raw_Reads/AOP39_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP39_D1_surf_bin.13.fa -o AOP39_D1_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP3_A1_surf_R1.fastq.gz /Raw_Reads/AOP3_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP3_A1_surf_bin.3.fa -o AOP3_A1_surf_bin.3
coverm genome --coupled /Raw_Reads/AOP3_F2_surf_R1.fastq.gz /Raw_Reads/AOP3_F2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP3_F2_surf_bin.1.fa -o AOP3_F2_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP40_10SA_surf_R1.fastq.gz /Raw_Reads/AOP40_10SA_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP40_10SA_surf_bin.23.fa -o AOP40_10SA_surf_bin.23
coverm genome --coupled /Raw_Reads/AOP40_10SA_surf_R1.fastq.gz /Raw_Reads/AOP40_10SA_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP40_10SA_surf_bin.2.fa -o AOP40_10SA_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP43_A2_surf_R1.fastq.gz /Raw_Reads/AOP43_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_A2_surf_bin.6.fa -o AOP43_A2_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP43_B1_surf_R1.fastq.gz /Raw_Reads/AOP43_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_B1_surf_bin.11.fa -o AOP43_B1_surf_bin.11
coverm genome --coupled /Raw_Reads/AOP43_B1_surf_R1.fastq.gz /Raw_Reads/AOP43_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_B1_surf_bin.6.fa -o AOP43_B1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP43_B2_surf_R1.fastq.gz /Raw_Reads/AOP43_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_B2_surf_bin.2.fa -o AOP43_B2_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP43_C2_coeur_R1.fastq.gz /Raw_Reads/AOP43_C2_coeur_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_C2_coeur_bin.16.fa -o AOP43_C2_coeur_bin.16
coverm genome --coupled /Raw_Reads/AOP43_C2_surf_R1.fastq.gz /Raw_Reads/AOP43_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_C2_surf_bin.11.fa -o AOP43_C2_surf_bin.11
coverm genome --coupled /Raw_Reads/AOP43_C2_surf_R1.fastq.gz /Raw_Reads/AOP43_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_C2_surf_bin.15.fa -o AOP43_C2_surf_bin.15
coverm genome --coupled /Raw_Reads/AOP43_C2_surf_R1.fastq.gz /Raw_Reads/AOP43_C2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_C2_surf_bin.26.fa -o AOP43_C2_surf_bin.26
coverm genome --coupled /Raw_Reads/AOP43_D1_surf_R1.fastq.gz /Raw_Reads/AOP43_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_D1_surf_bin.5.fa -o AOP43_D1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP43_E1_surf_R1.fastq.gz /Raw_Reads/AOP43_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_E1_surf_bin.10.fa -o AOP43_E1_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP43_E1_surf_R1.fastq.gz /Raw_Reads/AOP43_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_E1_surf_bin.12.fa -o AOP43_E1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP43_F2_surf_R1.fastq.gz /Raw_Reads/AOP43_F2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP43_F2_surf_bin.2.fa -o AOP43_F2_surf_bin.2
coverm genome --coupled /Raw_Reads/AOP44_A1_surf_R1.fastq.gz /Raw_Reads/AOP44_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP44_A1_surf_bin.6.fa -o AOP44_A1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP4_A1_surf_R1.fastq.gz /Raw_Reads/AOP4_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP4_A1_surf_bin.10.fa -o AOP4_A1_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP4_A1_surf_R1.fastq.gz /Raw_Reads/AOP4_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP4_A1_surf_bin.8.fa -o AOP4_A1_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP4_B2_surf_R1.fastq.gz /Raw_Reads/AOP4_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP4_B2_surf_bin.7.fa -o AOP4_B2_surf_bin.7
coverm genome --coupled /Raw_Reads/AOP4_B2_surf_R1.fastq.gz /Raw_Reads/AOP4_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP4_B2_surf_bin.8.fa -o AOP4_B2_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP4_E1_surf_R1.fastq.gz /Raw_Reads/AOP4_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP4_E1_surf_bin.13.fa -o AOP4_E1_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP4_E1_surf_R1.fastq.gz /Raw_Reads/AOP4_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP4_E1_surf_bin.1.fa -o AOP4_E1_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP5_A1_surf_R1.fastq.gz /Raw_Reads/AOP5_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_A1_surf_bin.19.fa -o AOP5_A1_surf_bin.19
coverm genome --coupled /Raw_Reads/AOP5_B1_surf_R1.fastq.gz /Raw_Reads/AOP5_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_B1_surf_bin.15.fa -o AOP5_B1_surf_bin.15
coverm genome --coupled /Raw_Reads/AOP5_B1_surf_R1.fastq.gz /Raw_Reads/AOP5_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_B1_surf_bin.24.fa -o AOP5_B1_surf_bin.24
coverm genome --coupled /Raw_Reads/AOP5_B1_surf_R1.fastq.gz /Raw_Reads/AOP5_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_B1_surf_bin.4.fa -o AOP5_B1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP5_B2_surf_R1.fastq.gz /Raw_Reads/AOP5_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_B2_surf_bin.3.fa -o AOP5_B2_surf_bin.3
coverm genome --coupled /Raw_Reads/AOP5_CZ1_surf_R1.fastq.gz /Raw_Reads/AOP5_CZ1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_CZ1_surf_bin.10.fa -o AOP5_CZ1_surf_bin.10
coverm genome --coupled /Raw_Reads/AOP5_CZ1_surf_R1.fastq.gz /Raw_Reads/AOP5_CZ1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_CZ1_surf_bin.12.fa -o AOP5_CZ1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP5_CZ1_surf_R1.fastq.gz /Raw_Reads/AOP5_CZ1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_CZ1_surf_bin.16.fa -o AOP5_CZ1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP5_FY1_surf_R1.fastq.gz /Raw_Reads/AOP5_FY1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_FY1_surf_bin.12.fa -o AOP5_FY1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP5_FY1_surf_R1.fastq.gz /Raw_Reads/AOP5_FY1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_FY1_surf_bin.4.fa -o AOP5_FY1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP5_FY1_surf_R1.fastq.gz /Raw_Reads/AOP5_FY1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_FY1_surf_bin.8.fa -o AOP5_FY1_surf_bin.8
coverm genome --coupled /Raw_Reads/AOP5_GZ1_surf_R1.fastq.gz /Raw_Reads/AOP5_GZ1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_GZ1_surf_bin.13.fa -o AOP5_GZ1_surf_bin.13
coverm genome --coupled /Raw_Reads/AOP5_GZ1_surf_R1.fastq.gz /Raw_Reads/AOP5_GZ1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP5_GZ1_surf_bin.21.fa -o AOP5_GZ1_surf_bin.21
coverm genome --coupled /Raw_Reads/AOP6_C1_surf_R1.fastq.gz /Raw_Reads/AOP6_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP6_C1_surf_bin.15.fa -o AOP6_C1_surf_bin.15
coverm genome --coupled /Raw_Reads/AOP6_C1_surf_R1.fastq.gz /Raw_Reads/AOP6_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP6_C1_surf_bin.1.fa -o AOP6_C1_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP7_A1_surf_R1.fastq.gz /Raw_Reads/AOP7_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_A1_surf_bin.12.fa -o AOP7_A1_surf_bin.12
coverm genome --coupled /Raw_Reads/AOP7_A1_surf_R1.fastq.gz /Raw_Reads/AOP7_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_A1_surf_bin.16.fa -o AOP7_A1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP7_A2_surf_R1.fastq.gz /Raw_Reads/AOP7_A2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_A2_surf_bin.9.fa -o AOP7_A2_surf_bin.9
coverm genome --coupled /Raw_Reads/AOP7_B1_surf_R1.fastq.gz /Raw_Reads/AOP7_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_B1_surf_bin.1.fa -o AOP7_B1_surf_bin.1
coverm genome --coupled /Raw_Reads/AOP7_B1_surf_R1.fastq.gz /Raw_Reads/AOP7_B1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_B1_surf_bin.3.fa -o AOP7_B1_surf_bin.3
coverm genome --coupled /Raw_Reads/AOP7_C1_surf_R1.fastq.gz /Raw_Reads/AOP7_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_C1_surf_bin.14.fa -o AOP7_C1_surf_bin.14
coverm genome --coupled /Raw_Reads/AOP7_C1_surf_R1.fastq.gz /Raw_Reads/AOP7_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_C1_surf_bin.16.fa -o AOP7_C1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP7_C1_surf_R1.fastq.gz /Raw_Reads/AOP7_C1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_C1_surf_bin.17.fa -o AOP7_C1_surf_bin.17
coverm genome --coupled /Raw_Reads/AOP7_D1_surf_R1.fastq.gz /Raw_Reads/AOP7_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_D1_surf_bin.16.fa -o AOP7_D1_surf_bin.16
coverm genome --coupled /Raw_Reads/AOP7_D1_surf_R1.fastq.gz /Raw_Reads/AOP7_D1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_D1_surf_bin.7.fa -o AOP7_D1_surf_bin.7
coverm genome --coupled /Raw_Reads/AOP7_E1_surf_R1.fastq.gz /Raw_Reads/AOP7_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_E1_surf_bin.11.fa -o AOP7_E1_surf_bin.11
coverm genome --coupled /Raw_Reads/AOP7_E1_surf_R1.fastq.gz /Raw_Reads/AOP7_E1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP7_E1_surf_bin.15.fa -o AOP7_E1_surf_bin.15
coverm genome --coupled /Raw_Reads/AOP8_A1_surf_R1.fastq.gz /Raw_Reads/AOP8_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP8_A1_surf_bin.11.fa -o AOP8_A1_surf_bin.11
coverm genome --coupled /Raw_Reads/AOP8_A1_surf_R1.fastq.gz /Raw_Reads/AOP8_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP8_A1_surf_bin.4.fa -o AOP8_A1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP8_B2_surf_R1.fastq.gz /Raw_Reads/AOP8_B2_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP8_B2_surf_bin.4.fa -o AOP8_B2_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP9_A1_surf_R1.fastq.gz /Raw_Reads/AOP9_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP9_A1_surf_bin.18.fa -o AOP9_A1_surf_bin.18
coverm genome --coupled /Raw_Reads/AOP9_A1_surf_R1.fastq.gz /Raw_Reads/AOP9_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP9_A1_surf_bin.4.fa -o AOP9_A1_surf_bin.4
coverm genome --coupled /Raw_Reads/AOP9_A1_surf_R1.fastq.gz /Raw_Reads/AOP9_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP9_A1_surf_bin.5.fa -o AOP9_A1_surf_bin.5
coverm genome --coupled /Raw_Reads/AOP9_A1_surf_R1.fastq.gz /Raw_Reads/AOP9_A1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP9_A1_surf_bin.6.fa -o AOP9_A1_surf_bin.6
coverm genome --coupled /Raw_Reads/AOP9_I1_surf_R1.fastq.gz /Raw_Reads/AOP9_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP9_I1_surf_bin.15.fa -o AOP9_I1_surf_bin.15
coverm genome --coupled /Raw_Reads/AOP9_I1_surf_R1.fastq.gz /Raw_Reads/AOP9_I1_surf_R2.fastq.gz --genome-fasta-files /Fungal_fasta/AOP9_I1_surf_bin.4.fa -o AOP9_I1_surf_bin.4

