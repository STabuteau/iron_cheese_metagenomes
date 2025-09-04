#Construct binary compressed datafiles for hmmscan and run hmmscan
conda activate hmmer-3.3.2
hmmpress All_HMMv4.hmm
for file in /METAGENOMES/*.faa;
do
   samplename=$(basename $file .faa)
   hmmscan --domtblout HMMER/${samplename}.tbl -E 0.00001 All_HMMv4.hmm /METAGENOMES/${samplename}.faa
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
done
rm -f *_parse_file.tmp
rm -f *_parse_merge.tsv

#Mapping of reads on hits to obtain quantitative data
conda activate seqkit-2.0.0
for file in /HMMER/*.tbl;
do
   samplename=$(basename $file .tbl)
   awk '{print $3}' ${samplename}_parse_merge_bitscore.tsv > ${samplename}_header_vf.txt
   seqkit grep -f ${samplename}_header_vf.txt /FASTA/${samplename}.ffn -o ${samplename}_hmmscan.fasta
   qsub -q short.q -wd /HMMER/ -pe thread 12 mapping_AOP.sh /RAW_READS_folder/ ${samplename}_hmmscan.fasta ${samplename}
done
rm -fr *_header_vf.txt

#Calculation of the weighted average length of hits recognised by each HMM
for file in /work_projet/metapdocheese/METAG/HMMSCAN/METAG_V3/*.tbl;
do
   samplename=$(basename $file .tbl)
   sed -i 's/#//' ${samplename}_cov.csv
   Rscript append_cov.R ${samplename}_parse_merge_bitscore.tsv ${samplename}_cov.csv ${samplename}_parse_merge_bitscore_cov.tsv
   sed -i '1i\HMMs_ID\thmm_len\tquery\tqlen\tbitscore\teval\thmm_from\thmm_to\tquery_from\tquery_to\tq_cov\tTarget\tIron acquisition system\tCategory\tSiderophore family\tImport componant\tHMMs_origin\tBitscores used in this study\tnumreads\tcoverage\tmeandepth' ${samplename}_parse_merge_bitscore_cov.tsv
   awk -F '\t' '{print $1"\t"$4"\t"$19}' ${samplename}_parse_merge_bitscore_cov.tsv > ${samplename}_parse_merge_bitscore_cov_crop.tsv
   awk 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}} END{ \
for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s" "a[j][i]}; print s}}' ${samplename}_parse_merge_bitscore_cov_crop.tsv > ${samplename}_parse_merge_bitscore_cov_crop_merged.tsv
   awk -F '\t' '{$4 = $2*$3 ; print}' ${samplename}_parse_merge_bitscore_cov_crop.tsv > ${samplename}_parse_merge_bitscore_cov_crop_1.tsv
   sed -i 's/ /\t/g' ${samplename}_parse_merge_bitscore_cov_crop_1.tsv
   awk -F '\t' 'NR==1{print;next} {for (i=2;i<=NF;i++) {a[$1][i]+=$i}} END{ \
for (j in a) {s=j; for (i=2;i<=NF;i++) {s=s" "a[j][i]}; print s}}' ${samplename}_parse_merge_bitscore_cov_crop_1.tsv > ${samplename}_parse_merge_bitscore_cov_crop_2.tsv
   sed -i -e "1d" ${samplename}_parse_merge_bitscore_cov_crop_2.tsv
   sed -i 's/ /\t/g' ${samplename}_parse_merge_bitscore_cov_crop_2.tsv
   awk '$3!=0' ${samplename}_parse_merge_bitscore_cov_crop_2.tsv | awk -F '\t' '{$5 = $4/$3 ; print}' | awk -F ' ' '{print $1 "\t" $5}'  > ${samplename}_parse_merge_bitscore_cov_crop_mean.tsv
   sed -i '1 i\HMMs_ID\tqlen_mean' ${samplename}_parse_merge_bitscore_cov_crop_mean.tsv
   sed -i 's/ /\t/g' ${samplename}_parse_merge_bitscore_cov_crop_mean.tsv
   sed -i 's/ /\t/g' ${samplename}_parse_merge_bitscore_cov_crop_merged.tsv
   Rscript merge_cov_qlen.R ${samplename}_parse_merge_bitscore_cov_crop_mean.tsv ${samplename}_parse_merge_bitscore_cov_crop_merged.tsv ${samplename}_parse_merge_bitscore_cov_crop_mean_merged.tsv
   sed '1d' ${samplename}_parse_merge_bitscore_cov_crop_mean_merged.tsv > ${samplename}_vf.tsv
done
rm -fr *_parse_merge_bitscore_cov_crop.tsv
rm -fr *_parse_merge_bitscore_cov_crop_merged.tsv
rm -fr *_parse_merge_bitscore_cov_crop_1.tsv
rm -fr *_parse_merge_bitscore_cov_crop_2.tsv
rm -fr *_parse_merge_bitscore_cov_crop_mean.tsv
rm -fr *_parse_merge_bitscore_cov_crop_mean_merged.tsv

#Normalization in RPKM (Reads Per Kilobase Million)
awk -F '\t' '{$4 = ($3/(17922643/1000000))/($2/1000) ; print}' AOP10_L1_surf_vf.tsv > AOP10_L1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(6765848/1000000))/($2/1000) ; print}' AOP10_LA_L1_vf.tsv > AOP10_LA_L1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19682603/1000000))/($2/1000) ; print}' AOP11_C1_surf_vf.tsv > AOP11_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(10533948/1000000))/($2/1000) ; print}' AOP12_A1_coeur_vf.tsv > AOP12_A1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18536302/1000000))/($2/1000) ; print}' AOP12_E2_surf_vf.tsv > AOP12_E2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(9481811/1000000))/($2/1000) ; print}' AOP13_A1_surf_vf.tsv > AOP13_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14829671/1000000))/($2/1000) ; print}' AOP13_A3_surf_vf.tsv > AOP13_A3_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14660971/1000000))/($2/1000) ; print}' AOP14_B2_coeur_vf.tsv > AOP14_B2_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14919336/1000000))/($2/1000) ; print}' AOP14_B2_surf_vf.tsv > AOP14_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14203540/1000000))/($2/1000) ; print}' AOP14_E2_surf_vf.tsv > AOP14_E2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13802250/1000000))/($2/1000) ; print}' AOP15_B1_surf_vf.tsv > AOP15_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12213855/1000000))/($2/1000) ; print}' AOP15_C2_surf_vf.tsv > AOP15_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16056017/1000000))/($2/1000) ; print}' AOP15_D2_surf_vf.tsv > AOP15_D2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17408302/1000000))/($2/1000) ; print}' AOP16_AQ2_surf_vf.tsv > AOP16_AQ2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17599572/1000000))/($2/1000) ; print}' AOP17_E1_surf_vf.tsv > AOP17_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15526321/1000000))/($2/1000) ; print}' AOP18_A1_surf_vf.tsv > AOP18_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13516254/1000000))/($2/1000) ; print}' AOP18_A2_surf_vf.tsv > AOP18_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15152193/1000000))/($2/1000) ; print}' AOP18_B2_surf_vf.tsv > AOP18_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(21062738/1000000))/($2/1000) ; print}' AOP18_D1_surf_vf.tsv > AOP18_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14364807/1000000))/($2/1000) ; print}' AOP18_LA_D1_vf.tsv > AOP18_LA_D1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13016677/1000000))/($2/1000) ; print}' AOP19_D1_surf_vf.tsv > AOP19_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17042297/1000000))/($2/1000) ; print}' AOP19_E1_surf_vf.tsv > AOP19_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17891608/1000000))/($2/1000) ; print}' AOP19_H1_surf_vf.tsv > AOP19_H1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17903799/1000000))/($2/1000) ; print}' AOP1_A1_surf_vf.tsv > AOP1_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16534874/1000000))/($2/1000) ; print}' AOP1_A2_coeur_vf.tsv > AOP1_A2_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14705038/1000000))/($2/1000) ; print}' AOP1_B1_coeur_vf.tsv > AOP1_B1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(20880597/1000000))/($2/1000) ; print}' AOP1_B1_surf_vf.tsv > AOP1_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(20231889/1000000))/($2/1000) ; print}' AOP1_C1_surf_vf.tsv > AOP1_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14297250/1000000))/($2/1000) ; print}' AOP1_D1_coeur_vf.tsv > AOP1_D1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17080611/1000000))/($2/1000) ; print}' AOP1_E1_surf_vf.tsv > AOP1_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12223988/1000000))/($2/1000) ; print}' AOP1_F1_coeur_vf.tsv > AOP1_F1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14303178/1000000))/($2/1000) ; print}' AOP1_F1_surf_vf.tsv > AOP1_F1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19195910/1000000))/($2/1000) ; print}' AOP1_G1_surf_vf.tsv > AOP1_G1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18309864/1000000))/($2/1000) ; print}' AOP1_H1_surf_vf.tsv > AOP1_H1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13241217/1000000))/($2/1000) ; print}' AOP1_I1_surf_vf.tsv > AOP1_I1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(4191243 /1000000))/($2/1000) ; print}' AOP1_LA_G1_vf.tsv > AOP1_LA_G1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18472474/1000000))/($2/1000) ; print}' AOP20_C1_surf_vf.tsv > AOP20_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15673535/1000000))/($2/1000) ; print}' AOP21_A1_surf_vf.tsv > AOP21_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16332973/1000000))/($2/1000) ; print}' AOP21_I1_surf_vf.tsv > AOP21_I1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(975112/1000000))/($2/1000) ; print}' AOP21_LA_D1_vf.tsv > AOP21_LA_D1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(11322202/1000000))/($2/1000) ; print}' AOP22_C1_surf_vf.tsv > AOP22_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(11795811/1000000))/($2/1000) ; print}' AOP22_E1_surf_vf.tsv > AOP22_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(2283380/1000000))/($2/1000) ; print}' AOP22_LA_C1_vf.tsv > AOP22_LA_C1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13579589/1000000))/($2/1000) ; print}' AOP23_AC1_coeur_vf.tsv > AOP23_AC1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16162653/1000000))/($2/1000) ; print}' AOP23_C1_coeur_vf.tsv > AOP23_C1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15588748/1000000))/($2/1000) ; print}' AOP23_C1_surf_vf.tsv > AOP23_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12410976/1000000))/($2/1000) ; print}' AOP23_I1_surf_vf.tsv > AOP23_I1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(11273706/1000000))/($2/1000) ; print}' AOP24_B1_coeur_vf.tsv > AOP24_B1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(8048925/1000000))/($2/1000) ; print}' AOP24_C2_surf_vf.tsv > AOP24_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(11872009/1000000))/($2/1000) ; print}' AOP24_D1_surf_vf.tsv > AOP24_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(3407977/1000000))/($2/1000) ; print}' AOP24_LA_D1_vf.tsv > AOP24_LA_D1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19317054/1000000))/($2/1000) ; print}' AOP25_A1_surf_vf.tsv > AOP25_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15683980/1000000))/($2/1000) ; print}' AOP25_B1_coeur_vf.tsv > AOP25_B1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15399951/1000000))/($2/1000) ; print}' AOP25_B1_surf_vf.tsv > AOP25_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13407781/1000000))/($2/1000) ; print}' AOP25_B2_surf_vf.tsv > AOP25_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18598046/1000000))/($2/1000) ; print}' AOP25_C1_surf_vf.tsv > AOP25_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18130009/1000000))/($2/1000) ; print}' AOP25_D2_surf_vf.tsv > AOP25_D2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(7342166/1000000))/($2/1000) ; print}' AOP25_F1_surf_vf.tsv > AOP25_F1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15057622/1000000))/($2/1000) ; print}' AOP25_F2_surf_vf.tsv > AOP25_F2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19237693/1000000))/($2/1000) ; print}' AOP26_D1_surf_vf.tsv > AOP26_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(10936158/1000000))/($2/1000) ; print}' AOP26_G1_coeur_vf.tsv > AOP26_G1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14598597/1000000))/($2/1000) ; print}' AOP27_A1_surf_vf.tsv > AOP27_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16203507/1000000))/($2/1000) ; print}' AOP28_A2_surf_vf.tsv > AOP28_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18169703/1000000))/($2/1000) ; print}' AOP28_C2_surf_vf.tsv > AOP28_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18486267/1000000))/($2/1000) ; print}' AOP28_E1_surf_vf.tsv > AOP28_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(9416593/1000000))/($2/1000) ; print}' AOP29_A1_coeur_vf.tsv > AOP29_A1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(27378003/1000000))/($2/1000) ; print}' AOP29_B2_surf_vf.tsv > AOP29_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(25397507/1000000))/($2/1000) ; print}' AOP29_E1_surf_vf.tsv > AOP29_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(10871746/1000000))/($2/1000) ; print}' AOP2_B2_coeur_vf.tsv > AOP2_B2_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(9764371 /1000000))/($2/1000) ; print}' AOP2_E1_surf_vf.tsv > AOP2_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15471367/1000000))/($2/1000) ; print}' AOP30_B1_surf_vf.tsv > AOP30_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19298220/1000000))/($2/1000) ; print}' AOP30_B2_surf_vf.tsv > AOP30_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17304720/1000000))/($2/1000) ; print}' AOP30_C1_coeur_vf.tsv > AOP30_C1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19277778/1000000))/($2/1000) ; print}' AOP30_C2_surf_vf.tsv > AOP30_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(21080649/1000000))/($2/1000) ; print}' AOP30_D2_surf_vf.tsv > AOP30_D2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(9555135/1000000))/($2/1000) ; print}' AOP31_A2_surf_vf.tsv > AOP31_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18656568/1000000))/($2/1000) ; print}' AOP31_B2_surf_vf.tsv > AOP31_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(20718464/1000000))/($2/1000) ; print}' AOP31_C2_surf_vf.tsv > AOP31_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(1888602 /1000000))/($2/1000) ; print}' AOP31_LA_C2_vf.tsv > AOP31_LA_C2_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17450595/1000000))/($2/1000) ; print}' AOP32_A1_surf_vf.tsv > AOP32_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13839319/1000000))/($2/1000) ; print}' AOP33_10_coeur_vf.tsv > AOP33_10_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17042130/1000000))/($2/1000) ; print}' AOP33_2_surf_vf.tsv > AOP33_2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(11305841/1000000))/($2/1000) ; print}' AOP33_4_coeur_vf.tsv > AOP33_4_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15712003/1000000))/($2/1000) ; print}' AOP33_8_surf_vf.tsv > AOP33_8_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(22175245/1000000))/($2/1000) ; print}' AOP34_AQ2_surf_vf.tsv > AOP34_AQ2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19094997/1000000))/($2/1000) ; print}' AOP34_BR2_surf_vf.tsv > AOP34_BR2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19930789/1000000))/($2/1000) ; print}' AOP34_CS2_surf_vf.tsv > AOP34_CS2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(21182312/1000000))/($2/1000) ; print}' AOP35_01H_surf_vf.tsv > AOP35_01H_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(21735110/1000000))/($2/1000) ; print}' AOP35_02H_surf_vf.tsv > AOP35_02H_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19232828/1000000))/($2/1000) ; print}' AOP35_03-ete_surf_vf.tsv > AOP35_03-ete_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(23286007/1000000))/($2/1000) ; print}' AOP35_03H_surf_vf.tsv > AOP35_03H_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(20090112/1000000))/($2/1000) ; print}' AOP35_04H_surf_vf.tsv > AOP35_04H_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(25274193/1000000))/($2/1000) ; print}' AOP35_05H_surf_vf.tsv > AOP35_05H_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17135926/1000000))/($2/1000) ; print}' AOP35_LA_E1_vf.tsv > AOP35_LA_E1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17055387/1000000))/($2/1000) ; print}' AOP36_A1_surf_vf.tsv > AOP36_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(24670242/1000000))/($2/1000) ; print}' AOP36_B1_surf_vf.tsv > AOP36_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16970914/1000000))/($2/1000) ; print}' AOP36_E1_surf_vf.tsv > AOP36_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14902780/1000000))/($2/1000) ; print}' AOP37_A2_surf_vf.tsv > AOP37_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18744912/1000000))/($2/1000) ; print}' AOP38_B2_surf_vf.tsv > AOP38_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13172983/1000000))/($2/1000) ; print}' AOP38_C1_surf_vf.tsv > AOP38_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15947066/1000000))/($2/1000) ; print}' AOP38_E1_surf_vf.tsv > AOP38_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(8845584/1000000))/($2/1000) ; print}' AOP39_D1_surf_vf.tsv > AOP39_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14537572/1000000))/($2/1000) ; print}' AOP3_A1_surf_vf.tsv > AOP3_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16079151/1000000))/($2/1000) ; print}' AOP3_F2_surf_vf.tsv > AOP3_F2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(20824379/1000000))/($2/1000) ; print}' AOP40_10SA_surf_vf.tsv > AOP40_10SA_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(10826731/1000000))/($2/1000) ; print}' AOP40_8SA_coeur_vf.tsv > AOP40_8SA_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13813859/1000000))/($2/1000) ; print}' AOP41_A1_coeur_vf.tsv > AOP41_A1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15484779/1000000))/($2/1000) ; print}' AOP42_A1_surf_vf.tsv > AOP42_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17017176/1000000))/($2/1000) ; print}' AOP42_A2_surf_vf.tsv > AOP42_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19380514/1000000))/($2/1000) ; print}' AOP42_C1_surf_vf.tsv > AOP42_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(22851995/1000000))/($2/1000) ; print}' AOP42_C2_surf_vf.tsv > AOP42_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15511476/1000000))/($2/1000) ; print}' AOP42_E1_surf_vf.tsv > AOP42_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17124377/1000000))/($2/1000) ; print}' AOP43_A2_surf_vf.tsv > AOP43_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15311235/1000000))/($2/1000) ; print}' AOP43_B1_surf_vf.tsv > AOP43_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(22806654/1000000))/($2/1000) ; print}' AOP43_B2_surf_vf.tsv > AOP43_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17948272/1000000))/($2/1000) ; print}' AOP43_C2_coeur_vf.tsv > AOP43_C2_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(23790322/1000000))/($2/1000) ; print}' AOP43_C2_surf_vf.tsv > AOP43_C2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(21353150/1000000))/($2/1000) ; print}' AOP43_D1_surf_vf.tsv > AOP43_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(20994617/1000000))/($2/1000) ; print}' AOP43_E1_surf_vf.tsv > AOP43_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(23397302/1000000))/($2/1000) ; print}' AOP43_F2_surf_vf.tsv > AOP43_F2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14644578/1000000))/($2/1000) ; print}' AOP44_A1_surf_vf.tsv > AOP44_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(1925080/1000000))/($2/1000) ; print}' AOP44_E1_coeur_vf.tsv > AOP44_E1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17067997/1000000))/($2/1000) ; print}' AOP4_A1_surf_vf.tsv > AOP4_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(6099106/1000000))/($2/1000) ; print}' AOP4_B2_surf_vf.tsv > AOP4_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17129472/1000000))/($2/1000) ; print}' AOP4_E1_surf_vf.tsv > AOP4_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12152825/1000000))/($2/1000) ; print}' AOP5_A1_surf_vf.tsv > AOP5_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19022872/1000000))/($2/1000) ; print}' AOP5_B1_surf_vf.tsv > AOP5_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(11959689/1000000))/($2/1000) ; print}' AOP5_B2_surf_vf.tsv > AOP5_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13174479/1000000))/($2/1000) ; print}' AOP5_CZ1_surf_vf.tsv > AOP5_CZ1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(4249377/1000000))/($2/1000) ; print}' AOP5_DX2_coeur_vf.tsv > AOP5_DX2_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(19155067/1000000))/($2/1000) ; print}' AOP5_FY1_surf_vf.tsv > AOP5_FY1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12273943/1000000))/($2/1000) ; print}' AOP5_GZ1_surf_vf.tsv > AOP5_GZ1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(8416420/1000000))/($2/1000) ; print}' AOP5_LA_FY1_vf.tsv > AOP5_LA_FY1_RPKM.tsv
awk -F '\t' '{$4 = ($3/(13239190/1000000))/($2/1000) ; print}' AOP6_C1_surf_vf.tsv > AOP6_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(17306982/1000000))/($2/1000) ; print}' AOP7_A1_surf_vf.tsv > AOP7_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(16988674/1000000))/($2/1000) ; print}' AOP7_A2_surf_vf.tsv > AOP7_A2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12456929/1000000))/($2/1000) ; print}' AOP7_B1_surf_vf.tsv > AOP7_B1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14047723/1000000))/($2/1000) ; print}' AOP7_C1_surf_vf.tsv > AOP7_C1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14015988/1000000))/($2/1000) ; print}' AOP7_D1_surf_vf.tsv > AOP7_D1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12646823/1000000))/($2/1000) ; print}' AOP7_D2_surf_vf.tsv > AOP7_D2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(12985374/1000000))/($2/1000) ; print}' AOP7_E1_surf_vf.tsv > AOP7_E1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14437611/1000000))/($2/1000) ; print}' AOP8_A1_surf_vf.tsv > AOP8_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(18917956/1000000))/($2/1000) ; print}' AOP8_B2_surf_vf.tsv > AOP8_B2_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14358394/1000000))/($2/1000) ; print}' AOP9_A1_surf_vf.tsv > AOP9_A1_surf_RPKM.tsv
awk -F '\t' '{$4 = ($3/(15914530/1000000))/($2/1000) ; print}' AOP9_B1_coeur_vf.tsv > AOP9_B1_coeur_RPKM.tsv
awk -F '\t' '{$4 = ($3/(14505529/1000000))/($2/1000) ; print}' AOP9_I1_surf_vf.tsv > AOP9_I1_surf_RPKM.tsv

#Preparation of a file with all data
for file in /HMMER/*.tbl;
do
   samplename=$(basename $file .tbl)
   sed -i 's/ /\t/g' ${samplename}_RPKM.tsv
   awk '{print $1 "\t" $4}' ${samplename}_RPKM.tsv > ${samplename}_RPKM_only.tsv
done
sed 's/$/\tAOP10_L1_surf/' AOP10_L1_surf_RPKM_only.tsv > AOP10_L1_surf_RPKM_vf.tsv
sed 's/$/\tAOP10_LA_L1/' AOP10_LA_L1_RPKM_only.tsv > AOP10_LA_L1_RPKM_vf.tsv
sed 's/$/\tAOP11_C1_surf/' AOP11_C1_surf_RPKM_only.tsv > AOP11_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP12_A1_coeur/' AOP12_A1_coeur_RPKM_only.tsv > AOP12_A1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP12_E2_surf/' AOP12_E2_surf_RPKM_only.tsv > AOP12_E2_surf_RPKM_vf.tsv
sed 's/$/\tAOP13_A1_surf/' AOP13_A1_surf_RPKM_only.tsv > AOP13_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP13_A3_surf/' AOP13_A3_surf_RPKM_only.tsv > AOP13_A3_surf_RPKM_vf.tsv
sed 's/$/\tAOP14_B2_coeur/' AOP14_B2_coeur_RPKM_only.tsv > AOP14_B2_coeur_RPKM_vf.tsv
sed 's/$/\tAOP14_B2_surf/' AOP14_B2_surf_RPKM_only.tsv > AOP14_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP14_E2_surf/' AOP14_E2_surf_RPKM_only.tsv > AOP14_E2_surf_RPKM_vf.tsv
sed 's/$/\tAOP15_B1_surf/' AOP15_B1_surf_RPKM_only.tsv > AOP15_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP15_C2_surf/' AOP15_C2_surf_RPKM_only.tsv > AOP15_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP15_D2_surf/' AOP15_D2_surf_RPKM_only.tsv > AOP15_D2_surf_RPKM_vf.tsv
sed 's/$/\tAOP16_AQ2_surf/' AOP16_AQ2_surf_RPKM_only.tsv > AOP16_AQ2_surf_RPKM_vf.tsv
sed 's/$/\tAOP17_E1_surf/' AOP17_E1_surf_RPKM_only.tsv > AOP17_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP18_A1_surf/' AOP18_A1_surf_RPKM_only.tsv > AOP18_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP18_A2_surf/' AOP18_A2_surf_RPKM_only.tsv > AOP18_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP18_B2_surf/' AOP18_B2_surf_RPKM_only.tsv > AOP18_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP18_D1_surf/' AOP18_D1_surf_RPKM_only.tsv > AOP18_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP18_LA_D1/' AOP18_LA_D1_RPKM_only.tsv > AOP18_LA_D1_RPKM_vf.tsv
sed 's/$/\tAOP19_D1_surf/' AOP19_D1_surf_RPKM_only.tsv > AOP19_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP19_E1_surf/' AOP19_E1_surf_RPKM_only.tsv > AOP19_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP19_H1_surf/' AOP19_H1_surf_RPKM_only.tsv > AOP19_H1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_A1_surf/' AOP1_A1_surf_RPKM_only.tsv > AOP1_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_A2_coeur/' AOP1_A2_coeur_RPKM_only.tsv > AOP1_A2_coeur_RPKM_vf.tsv
sed 's/$/\tAOP1_B1_coeur/' AOP1_B1_coeur_RPKM_only.tsv > AOP1_B1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP1_B1_surf/' AOP1_B1_surf_RPKM_only.tsv > AOP1_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_C1_surf/' AOP1_C1_surf_RPKM_only.tsv > AOP1_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_D1_coeur/' AOP1_D1_coeur_RPKM_only.tsv > AOP1_D1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP1_E1_surf/' AOP1_E1_surf_RPKM_only.tsv > AOP1_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_F1_coeur/' AOP1_F1_coeur_RPKM_only.tsv > AOP1_F1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP1_F1_surf/' AOP1_F1_surf_RPKM_only.tsv > AOP1_F1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_G1_surf/' AOP1_G1_surf_RPKM_only.tsv > AOP1_G1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_H1_surf/' AOP1_H1_surf_RPKM_only.tsv > AOP1_H1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_I1_surf/' AOP1_I1_surf_RPKM_only.tsv > AOP1_I1_surf_RPKM_vf.tsv
sed 's/$/\tAOP1_LA_G1/' AOP1_LA_G1_RPKM_only.tsv > AOP1_LA_G1_RPKM_vf.tsv
sed 's/$/\tAOP20_C1_surf/' AOP20_C1_surf_RPKM_only.tsv > AOP20_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP21_A1_surf/' AOP21_A1_surf_RPKM_only.tsv > AOP21_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP21_I1_surf/' AOP21_I1_surf_RPKM_only.tsv > AOP21_I1_surf_RPKM_vf.tsv
sed 's/$/\tAOP21_LA_D1/' AOP21_LA_D1_RPKM_only.tsv > AOP21_LA_D1_RPKM_vf.tsv
sed 's/$/\tAOP22_C1_surf/' AOP22_C1_surf_RPKM_only.tsv > AOP22_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP22_E1_surf/' AOP22_E1_surf_RPKM_only.tsv > AOP22_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP22_LA_C1/' AOP22_LA_C1_RPKM_only.tsv > AOP22_LA_C1_RPKM_vf.tsv
sed 's/$/\tAOP23_AC1_coeur/' AOP23_AC1_coeur_RPKM_only.tsv > AOP23_AC1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP23_C1_coeur/' AOP23_C1_coeur_RPKM_only.tsv > AOP23_C1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP23_C1_surf/' AOP23_C1_surf_RPKM_only.tsv > AOP23_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP23_I1_surf/' AOP23_I1_surf_RPKM_only.tsv > AOP23_I1_surf_RPKM_vf.tsv
sed 's/$/\tAOP24_B1_coeur/' AOP24_B1_coeur_RPKM_only.tsv > AOP24_B1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP24_C2_surf/' AOP24_C2_surf_RPKM_only.tsv > AOP24_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP24_D1_surf/' AOP24_D1_surf_RPKM_only.tsv > AOP24_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP24_LA_D1/' AOP24_LA_D1_RPKM_only.tsv > AOP24_LA_D1_RPKM_vf.tsv
sed 's/$/\tAOP25_A1_surf/' AOP25_A1_surf_RPKM_only.tsv > AOP25_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP25_B1_coeur/' AOP25_B1_coeur_RPKM_only.tsv > AOP25_B1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP25_B1_surf/' AOP25_B1_surf_RPKM_only.tsv > AOP25_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP25_B2_surf/' AOP25_B2_surf_RPKM_only.tsv > AOP25_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP25_C1_surf/' AOP25_C1_surf_RPKM_only.tsv > AOP25_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP25_D2_surf/' AOP25_D2_surf_RPKM_only.tsv > AOP25_D2_surf_RPKM_vf.tsv
sed 's/$/\tAOP25_F1_surf/' AOP25_F1_surf_RPKM_only.tsv > AOP25_F1_surf_RPKM_vf.tsv
sed 's/$/\tAOP25_F2_surf/' AOP25_F2_surf_RPKM_only.tsv > AOP25_F2_surf_RPKM_vf.tsv
sed 's/$/\tAOP26_D1_surf/' AOP26_D1_surf_RPKM_only.tsv > AOP26_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP26_G1_coeur/' AOP26_G1_coeur_RPKM_only.tsv > AOP26_G1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP27_A1_surf/' AOP27_A1_surf_RPKM_only.tsv > AOP27_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP28_A2_surf/' AOP28_A2_surf_RPKM_only.tsv > AOP28_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP28_C2_surf/' AOP28_C2_surf_RPKM_only.tsv > AOP28_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP28_E1_surf/' AOP28_E1_surf_RPKM_only.tsv > AOP28_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP29_A1_coeur/' AOP29_A1_coeur_RPKM_only.tsv > AOP29_A1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP29_B2_surf/' AOP29_B2_surf_RPKM_only.tsv > AOP29_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP29_E1_surf/' AOP29_E1_surf_RPKM_only.tsv > AOP29_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP2_B2_coeur/' AOP2_B2_coeur_RPKM_only.tsv > AOP2_B2_coeur_RPKM_vf.tsv
sed 's/$/\tAOP2_E1_surf/' AOP2_E1_surf_RPKM_only.tsv > AOP2_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP30_B1_surf/' AOP30_B1_surf_RPKM_only.tsv > AOP30_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP30_B2_surf/' AOP30_B2_surf_RPKM_only.tsv > AOP30_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP30_C1_coeur/' AOP30_C1_coeur_RPKM_only.tsv > AOP30_C1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP30_C2_surf/' AOP30_C2_surf_RPKM_only.tsv > AOP30_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP30_D2_surf/' AOP30_D2_surf_RPKM_only.tsv > AOP30_D2_surf_RPKM_vf.tsv
sed 's/$/\tAOP31_A2_surf/' AOP31_A2_surf_RPKM_only.tsv > AOP31_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP31_B2_surf/' AOP31_B2_surf_RPKM_only.tsv > AOP31_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP31_C2_surf/' AOP31_C2_surf_RPKM_only.tsv > AOP31_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP31_LA_C2/' AOP31_LA_C2_RPKM_only.tsv > AOP31_LA_C2_RPKM_vf.tsv
sed 's/$/\tAOP32_A1_surf/' AOP32_A1_surf_RPKM_only.tsv > AOP32_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP33_10_coeur/' AOP33_10_coeur_RPKM_only.tsv > AOP33_10_coeur_RPKM_vf.tsv
sed 's/$/\tAOP33_2_surf/' AOP33_2_surf_RPKM_only.tsv > AOP33_2_surf_RPKM_vf.tsv
sed 's/$/\tAOP33_4_coeur/' AOP33_4_coeur_RPKM_only.tsv > AOP33_4_coeur_RPKM_vf.tsv
sed 's/$/\tAOP33_8_surf/' AOP33_8_surf_RPKM_only.tsv > AOP33_8_surf_RPKM_vf.tsv
sed 's/$/\tAOP34_AQ2_surf/' AOP34_AQ2_surf_RPKM_only.tsv > AOP34_AQ2_surf_RPKM_vf.tsv
sed 's/$/\tAOP34_BR2_surf/' AOP34_BR2_surf_RPKM_only.tsv > AOP34_BR2_surf_RPKM_vf.tsv
sed 's/$/\tAOP34_CS2_surf/' AOP34_CS2_surf_RPKM_only.tsv > AOP34_CS2_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_01H_surf/' AOP35_01H_surf_RPKM_only.tsv > AOP35_01H_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_02H_surf/' AOP35_02H_surf_RPKM_only.tsv > AOP35_02H_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_03-ete_surf/' AOP35_03-ete_surf_RPKM_only.tsv > AOP35_03-ete_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_03H_surf/' AOP35_03H_surf_RPKM_only.tsv > AOP35_03H_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_04H_surf/' AOP35_04H_surf_RPKM_only.tsv > AOP35_04H_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_05H_surf/' AOP35_05H_surf_RPKM_only.tsv > AOP35_05H_surf_RPKM_vf.tsv
sed 's/$/\tAOP35_LA_E1/' AOP35_LA_E1_RPKM_only.tsv > AOP35_LA_E1_RPKM_vf.tsv
sed 's/$/\tAOP36_A1_surf/' AOP36_A1_surf_RPKM_only.tsv > AOP36_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP36_B1_surf/' AOP36_B1_surf_RPKM_only.tsv > AOP36_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP36_E1_surf/' AOP36_E1_surf_RPKM_only.tsv > AOP36_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP37_A2_surf/' AOP37_A2_surf_RPKM_only.tsv > AOP37_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP38_B2_surf/' AOP38_B2_surf_RPKM_only.tsv > AOP38_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP38_C1_surf/' AOP38_C1_surf_RPKM_only.tsv > AOP38_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP38_E1_surf/' AOP38_E1_surf_RPKM_only.tsv > AOP38_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP39_D1_surf/' AOP39_D1_surf_RPKM_only.tsv > AOP39_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP3_A1_surf/' AOP3_A1_surf_RPKM_only.tsv > AOP3_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP3_F2_surf/' AOP3_F2_surf_RPKM_only.tsv > AOP3_F2_surf_RPKM_vf.tsv
sed 's/$/\tAOP40_10SA_surf/' AOP40_10SA_surf_RPKM_only.tsv > AOP40_10SA_surf_RPKM_vf.tsv
sed 's/$/\tAOP40_8SA_coeur/' AOP40_8SA_coeur_RPKM_only.tsv > AOP40_8SA_coeur_RPKM_vf.tsv
sed 's/$/\tAOP41_A1_coeur/' AOP41_A1_coeur_RPKM_only.tsv > AOP41_A1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP42_A1_surf/' AOP42_A1_surf_RPKM_only.tsv > AOP42_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP42_A2_surf/' AOP42_A2_surf_RPKM_only.tsv > AOP42_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP42_C1_surf/' AOP42_C1_surf_RPKM_only.tsv > AOP42_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP42_C2_surf/' AOP42_C2_surf_RPKM_only.tsv > AOP42_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP42_E1_surf/' AOP42_E1_surf_RPKM_only.tsv > AOP42_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_A2_surf/' AOP43_A2_surf_RPKM_only.tsv > AOP43_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_B1_surf/' AOP43_B1_surf_RPKM_only.tsv > AOP43_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_B2_surf/' AOP43_B2_surf_RPKM_only.tsv > AOP43_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_C2_coeur/' AOP43_C2_coeur_RPKM_only.tsv > AOP43_C2_coeur_RPKM_vf.tsv
sed 's/$/\tAOP43_C2_surf/' AOP43_C2_surf_RPKM_only.tsv > AOP43_C2_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_D1_surf/' AOP43_D1_surf_RPKM_only.tsv > AOP43_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_E1_surf/' AOP43_E1_surf_RPKM_only.tsv > AOP43_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP43_F2_surf/' AOP43_F2_surf_RPKM_only.tsv > AOP43_F2_surf_RPKM_vf.tsv
sed 's/$/\tAOP44_A1_surf/' AOP44_A1_surf_RPKM_only.tsv > AOP44_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP44_E1_coeur/' AOP44_E1_coeur_RPKM_only.tsv > AOP44_E1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP4_A1_surf/' AOP4_A1_surf_RPKM_only.tsv > AOP4_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP4_B2_surf/' AOP4_B2_surf_RPKM_only.tsv > AOP4_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP4_E1_surf/' AOP4_E1_surf_RPKM_only.tsv > AOP4_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_A1_surf/' AOP5_A1_surf_RPKM_only.tsv > AOP5_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_B1_surf/' AOP5_B1_surf_RPKM_only.tsv > AOP5_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_B2_surf/' AOP5_B2_surf_RPKM_only.tsv > AOP5_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_CZ1_surf/' AOP5_CZ1_surf_RPKM_only.tsv > AOP5_CZ1_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_DX2_coeur/' AOP5_DX2_coeur_RPKM_only.tsv > AOP5_DX2_coeur_RPKM_vf.tsv
sed 's/$/\tAOP5_FY1_surf/' AOP5_FY1_surf_RPKM_only.tsv > AOP5_FY1_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_GZ1_surf/' AOP5_GZ1_surf_RPKM_only.tsv > AOP5_GZ1_surf_RPKM_vf.tsv
sed 's/$/\tAOP5_LA_FY1/' AOP5_LA_FY1_RPKM_only.tsv > AOP5_LA_FY1_RPKM_vf.tsv
sed 's/$/\tAOP6_C1_surf/' AOP6_C1_surf_RPKM_only.tsv > AOP6_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_A1_surf/' AOP7_A1_surf_RPKM_only.tsv > AOP7_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_A2_surf/' AOP7_A2_surf_RPKM_only.tsv > AOP7_A2_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_B1_surf/' AOP7_B1_surf_RPKM_only.tsv > AOP7_B1_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_C1_surf/' AOP7_C1_surf_RPKM_only.tsv > AOP7_C1_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_D1_surf/' AOP7_D1_surf_RPKM_only.tsv > AOP7_D1_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_D2_surf/' AOP7_D2_surf_RPKM_only.tsv > AOP7_D2_surf_RPKM_vf.tsv
sed 's/$/\tAOP7_E1_surf/' AOP7_E1_surf_RPKM_only.tsv > AOP7_E1_surf_RPKM_vf.tsv
sed 's/$/\tAOP8_A1_surf/' AOP8_A1_surf_RPKM_only.tsv > AOP8_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP8_B2_surf/' AOP8_B2_surf_RPKM_only.tsv > AOP8_B2_surf_RPKM_vf.tsv
sed 's/$/\tAOP9_A1_surf/' AOP9_A1_surf_RPKM_only.tsv > AOP9_A1_surf_RPKM_vf.tsv
sed 's/$/\tAOP9_B1_coeur/' AOP9_B1_coeur_RPKM_only.tsv > AOP9_B1_coeur_RPKM_vf.tsv
sed 's/$/\tAOP9_I1_surf/' AOP9_I1_surf_RPKM_only.tsv > AOP9_I1_surf_RPKM_vf.tsv
cat *_RPKM_vf.tsv > merged_RPKM.tsv
sed -i '1 i\HMMs_ID\tRPKM\tMetagenomes' merged_RPKM.tsv