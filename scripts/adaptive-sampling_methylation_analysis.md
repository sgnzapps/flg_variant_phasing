# FLG adaptive sampling methylation analysis

## Get FLG gene feature coordinates from GTF file

```bash
cd /public-data/references/GRCh38
zcat gencode.v43.primary_assembly.basic.annotation.gtf.gz | grep 'FLG"'

"""
chr1:152325189-152325239    exon1
chr1:152315319-152315477    exon2
chr1:152302165-152314747    exon3
chr1:152315478-152325188    intron1
chr1:152314748-152315318    intron2
chr1:152325240-152326239    promotor
"""
```

## Run modkit

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis
mkdir modkit
cd modkit
ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152325240-152326239 --cpg -r ${ref} ../whatshap/S11_11.haplotagged.bam ./S11_11.FLG_promo.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152325189-152325239 --cpg -r ${ref} ../whatshap/S11_11.haplotagged.bam ./S11_11.FLG_exon1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152315478-152325188 --cpg -r ${ref} ../whatshap/S11_11.haplotagged.bam ./S11_11.FLG_intron1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152315319-152315477 --cpg -r ${ref} ../whatshap/S11_11.haplotagged.bam ./S11_11.FLG_exon2.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152314748-152315318 --cpg -r ${ref} ../whatshap/S11_11.haplotagged.bam ./S11_11.FLG_intron2.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152302165-152314747 --cpg -r ${ref} ../whatshap/S11_11.haplotagged.bam ./S11_11.FLG_exon3.bed

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis
mkdir modkit
cd modkit
ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152325240-152326239 --cpg -r ${ref} ../whatshap/S12_12.haplotagged.bam ./S12_12.FLG_promo.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152325189-152325239 --cpg -r ${ref} ../whatshap/S12_12.haplotagged.bam ./S12_12.FLG_exon1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152315478-152325188 --cpg -r ${ref} ../whatshap/S12_12.haplotagged.bam ./S12_12.FLG_intron1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152315319-152315477 --cpg -r ${ref} ../whatshap/S12_12.haplotagged.bam ./S12_12.FLG_exon2.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152314748-152315318 --cpg -r ${ref} ../whatshap/S12_12.haplotagged.bam ./S12_12.FLG_intron2.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region chr1:152302165-152314747 --cpg -r ${ref} ../whatshap/S12_12.haplotagged.bam ./S12_12.FLG_exon3.bed

```

## Calculate methylation percentage for gene features

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis/modkit
COV=4  # Set Nvalidcov threshold as min 4

touch output.csv
truncate -s 0 output.csv  # ensure file is empty
while read s; do 
    i="S11_11.FLG_"${s}
    echo ${i} | sed 's/.bed//g' | sed 's/S11_11.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
    echo ${i}"-" | sed 's/.bed//g' | sed 's/S11_11.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp 
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
done < order.txt
rm *.tmp

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis/modkit
COV=4

touch output.csv
truncate -s 0 output.csv  # ensure file is empty
while read s; do 
    i="S12_12.FLG_"${s}
    echo ${i} | sed 's/.bed//g' | sed 's/S12_12.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
    echo ${i}"-" | sed 's/.bed//g' | sed 's/S12_12.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp 
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
done < order.txt
rm *.tmp

```

## Map uBAM to 11/12 rep reference

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis

samtools view S11_11_pass_F2308_q60_FLG.bam -H > header.tmp
samtools view S11_11_pass_F2308_q60_FLG.bam > body.tmp
nano header.tmp  # Add SM tag in @RG and additional @RG 70ee34a18a86cd9b21d361e75b4cfea5e8f19407_2021-05-17_dna_r9.4.1_minion_768_2f1c8637
cat header.tmp body.tmp | samtools view -Sb - -o S11_11_pass_F2308_q60_FLG_sm.bam

# Convert BAM to uBAM
java -jar ~/programs/picard.jar RevertSam I=S11_11_pass_F2308_q60_FLG_sm.bam O=S11_11_pass_F2308_q60_FLG_sm_umap.bam

REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/analysis_1223/FLG_11-RPT-8.1p.fa'
dorado_0.5.0 aligner -t 24 ${REF} S11_11_pass_F2308_q60_FLG_sm_umap.bam | samtools sort - -o S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam
samtools index S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis

samtools view S12_12_pass_F2308_q60_FLG.bam -H > header.tmp
samtools view S12_12_pass_F2308_q60_FLG.bam > body.tmp
nano header.tmp  # Add SM tag in @RG and additional @RG a907b8b8d679724cd450d5f5b8a616990ee9700c_2021-05-17_dna_r9.4.1_minion_768_2f1c8637
cat header.tmp body.tmp | samtools view -Sb - -o S12_12_pass_F2308_q60_FLG_sm.bam

# Convert BAM to uBAM
java -jar ~/programs/picard.jar RevertSam I=S12_12_pass_F2308_q60_FLG_sm.bam O=S12_12_pass_F2308_q60_FLG_sm_umap.bam

REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_12-RPT-8.1p-10.1p.fa'
dorado_0.5.0 aligner -t 24 ${REF} S12_12_pass_F2308_q60_FLG_sm_umap.bam | samtools sort - -o S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam
samtools index S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam
```

## Get RPT coordinates from 10.1038/ng2020

```bash
"""
# Refer to FLG_RPT_coord_sandilands.xlsx (Inferred from Sandilands et al., 2007, sup fig 1)

# RPT coordinates for FLG_11-RPT-8.1p.fa
S-100 FLG_11-RPT-8.1p:15393-16428
RPT0 FLG_11-RPT-8.1p:14478-15392
RPT1 FLG_11-RPT-8.1p:13503-14477
RPT2 FLG_11-RPT-8.1p:12531-13502
RPT3 FLG_11-RPT-8.1p:11559-12530
RPT4 FLG_11-RPT-8.1p:10586-11558
RPT5 FLG_11-RPT-8.1p:9612-10585
RPT6 FLG_11-RPT-8.1p:8640-9611
RPT7 FLG_11-RPT-8.1p:7668-8639
RPT8.1 FLG_11-RPT-8.1p:6696-7667
RPT8.1p FLG_11-RPT-8.1p:5724-6695
RPT9 FLG_11-RPT-8.1p:4752-5723
RPT10.1 FLG_11-RPT-8.1p:3780-4751
RPT11 FLG_11-RPT-8.1p:2703-3779
3UTR FLG_11-RPT-8.1p:2165-2702

# RPT coordinates for FLG_12-RPT-8.1p-10.1p.fa
S-100 FLG_12-RPT-8.1p-10.1p:16365-17400
RPT0 FLG_12-RPT-8.1p-10.1p:15450-16364
RPT1 FLG_12-RPT-8.1p-10.1p:14475-15449
RPT2 FLG_12-RPT-8.1p-10.1p:13503-14474
RPT3 FLG_12-RPT-8.1p-10.1p:12531-13502
RPT4 FLG_12-RPT-8.1p-10.1p:11558-12530
RPT5 FLG_12-RPT-8.1p-10.1p:10584-11557
RPT6 FLG_12-RPT-8.1p-10.1p:9612-10583
RPT7 FLG_12-RPT-8.1p-10.1p:8640-9611
RPT8.1 FLG_12-RPT-8.1p-10.1p:7668-8639
RPT8.1p FLG_12-RPT-8.1p-10.1p:6696-7667
RPT9 FLG_12-RPT-8.1p-10.1p:5724-6695
RPT10.1 FLG_12-RPT-8.1p-10.1p:4752-5723
RPT10.1p FLG_12-RPT-8.1p-10.1p:3780-4751
RPT11 FLG_12-RPT-8.1p-10.1p:2703-3779
3UTR FLG_12-RPT-8.1p-10.1p:2165-2702
"""
```

## Run modkit for each RPT

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis
mkdir rpt_mods
REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/analysis_1223/FLG_11-RPT-8.1p.fa'
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:15393-16428 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_s100.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:14478-15392 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt0.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:13503-14477 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:12531-13502 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt2.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:11559-12530 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt3.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:10586-11558 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt4.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:9612-10585 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt5.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:8640-9611 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt6.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:7668-8639 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt7.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:6696-7667 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt8.1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:5724-6695 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt8.1p.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:4752-5723 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt9.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:3780-4751 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt10.1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:2703-3779 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_rpt11.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_11-RPT-8.1p:2165-2702 --cpg -r ${REF} ./S11_11_pass_F2308_q60_FLG_sm_11ref-mods.bam ./rpt_mods/S11_11.FLG_3utr.bed

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis
mkdir rpt_mods
REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_12-RPT-8.1p-10.1p.fa'
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:16365-17400 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_s100.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:15450-16364 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt0.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:14475-15449 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:13503-14474 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt2.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:12531-13502 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt3.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:11558-12530 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt4.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:10584-11557 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt5.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:9612-10583 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt6.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:8640-9611 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt7.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:7668-8639 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt8.1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:6696-7667 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt8.1p.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:5724-6695 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt9.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:4752-5723 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt10.1.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:3780-4751 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt10.1p.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:2703-3779 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_rpt11.bed
~/programs/modkit/modkit pileup -t 24 --filter-threshold 0.75 --region FLG_12-RPT-8.1p-10.1p:2165-2702 --cpg -r ${REF} ./S12_12_pass_F2308_q60_FLG_sm_12ref-mods.bam ./rpt_mods/S12_12.FLG_3utr.bed

```

## Calculate methylation percentage for each RPT

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis/rpt_mods
COV=4
touch output.csv
truncate -s 0 output.csv  # ensure file is empty
while read s; do 
    i="S11_11.FLG_"${s}
    echo ${i} | sed 's/.bed//g' | sed 's/S11_11.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
    echo ${i}"-" | sed 's/.bed//g' | sed 's/S11_11.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp 
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
done < order.txt
rm *.tmp

# Just for double checking methylation percentage
for i in *.bed; do
    echo ${i}
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1,4 -d ' ' | awk '{s+=$1}{j+=$2}END{print (1-(j/s))*100}'
    echo ""
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1,4 -d ' ' | awk '{s+=$1}{j+=$2}END{print (1-(j/s))*100}'
    echo ""
done

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis/rpt_mods
COV=4   
touch output.csv
truncate -s 0 output.csv  # ensure file is empty
while read s; do 
    i="S12_12.FLG_"${s}
    echo ${i} | sed 's/.bed//g' | sed 's/S12_12.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="+") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
    echo ${i}"-" | sed 's/.bed//g' | sed 's/S12_12.FLG_//g' > ${i}.tmp
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | awk 'END{print NR}' >> ${i}.tmp # CpG count
    echo "" >> ${i}.tmp 
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 1 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N valid cov
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5mC
    awk -v var="${COV}" '{if (($4=="h") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 3 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N 5hmC
    awk -v var="${COV}" '{if (($4=="m") && ($6=="-") && ($5>=var)) print $0}' ${i} | cut -f 10 | cut -f 4 -d ' ' | awk '{s+=$1}END{print s}' >> ${i}.tmp  # N canonical
    paste -d "," output.csv ${i}.tmp > output2.csv
    mv output2.csv output.csv
done < order.txt
rm *.tmp

# Just for double checking methylation percentage
for i in *.bed; do
    echo ${i}
    awk '{if (($4=="m") && ($6=="+")) print $0}' ${i} | cut -f 10 | cut -f 1,4 -d ' ' | awk '{s+=$1}{j+=$2}END{print (1-(j/s))*100}'
    echo ""
    awk '{if (($4=="m") && ($6=="-")) print $0}' ${i} | cut -f 10 | cut -f 1,4 -d ' ' | awk '{s+=$1}{j+=$2}END{print (1-(j/s))*100}'
    echo ""
done
```
