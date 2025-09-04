#Construct binary compressed datafiles for hmmscan and run hmmscan

conda activate hmmer-3.3.2
hmmpress All_HMMv4.hmm
for file in /MAGs_GENOMES/*.faa;
do
   samplename=$(basename $file .faa)
   hmmscan --domtblout HMMER/${samplename}.tbl -E 0.00001 All_HMMv4.hmm /MAGs_GENOMES/${samplename}.faa
done
conda deactivate

#Filtered by predefined bitscore cutoffs for each HMM and length of the CDS recognised by the HMMs
for file in /HMMER/*.tbl;
do
   samplename=$(basename $file .tbl)
   cat ${samplename}.tbl | grep -v '^#' | \
   sed 's/\s\s*/\t/g' | \
   sort -k4,4 | \
   awk 'BEGIN{FS="\t"; OFS="\t"}{if(($19-$18+1)/$6>=0.7){print $0,($19-$18+1),($19-$18+1)/$6}}' | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,$3,$4,$6,$8,$13,$16,$17,$18,$19,$NF}' | \
   sed 's/,/\./g' > ${samplename}_parse_file.tmp
   Rscript merge_files.R ${samplename}_parse_file7.tmp HMMs_metadata.tsv ${samplename}_parse_merge.tsv
   awk -F '\t' '{if($5>=$NF) {print $0}}' ${samplename}_parse_merge.tsv | awk '$5 > max[$3] { max[$3] = $5; m[$3] = $0 } END { for (i in m) { print m[i] } }' > ${samplename}_parse_merge_bitscore.tsv
   awk -F "\t" '{print $1}' ${samplename}_parse_merge_bitscore.tsv | sed 's/$/\t1/' | sort | sed '1iHMMs_ID\tCount' | awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}} END{ for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s" "a[j][i]}; print s}}' | sed 's/ /\t/' | sed 's/$/\t'${samplename}'/' | sed '1d' | sed '1iHMMs_ID\tCount\tGenomes' > ${samplename}_parse_merge_bitscore_Count.tsv
done
rm -f *_parse_file.tmp
rm -f *_parse_merge.tsv
cat *_parse_merge_bitscore_Count.tsv | grep -v 'Count' | sed '1iHMMs_ID\tCount\tGenomes' > Genomes_parse_merge_bitscore_Count.tsv