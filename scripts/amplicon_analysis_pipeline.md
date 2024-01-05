# Filaggrin amplicon variant calling and phasing bioinformatics pipeline

## Pipeline

### Sample ID legend

```python
sample_legend = {
    '1A': 'R08_b01', 
    '2A': 'R08_b02', 
    '3A': 'R08_b03', 
    '4A': 'R08_b04', 
    '5A': 'R08_b05', 
    '6A': 'R08_b06', 
    '1B': 'R37_b01', 
    '2B': 'R37_b04', 
    '3B': 'R37_b05', 
    '4B',: 'R37_b07', 
    '5B': 'R37_b08', 
    '6B': 'R37_b11', 
    '7B': 'R37_b12', 
    '8B': 'R37_b02', 
    '9B': 'R37_b03', 
    '10B': 'R37_b06', 
    '11B': 'R37_b09', 
    '12B': 'R37_b10', 
    '13B': 'R37_b13', 
    '14B': 'R37_b14', 
    '15B': 'R37_b15', 
    '16B': 'R37_b16'
}
```

### Concatenate all passed FASTQ reads

```bash
cd /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5
mkdir analysis_1223 && cd analysis_1223
mkdir fastq && cd fastq
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00008_fast5/pass/barcode01/*.fastq.gz > R08_b01.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00008_fast5/pass/barcode02/*.fastq.gz > R08_b02.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00008_fast5/pass/barcode03/*.fastq.gz > R08_b03.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00008_fast5/pass/barcode04/*.fastq.gz > R08_b04.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00008_fast5/pass/barcode05/*.fastq.gz > R08_b05.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00008_fast5/pass/barcode06/*.fastq.gz > R08_b06.fq.gz

cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode01/*.fastq.gz > R37_b01.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode02/*.fastq.gz > R37_b02.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode03/*.fastq.gz > R37_b03.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode04/*.fastq.gz > R37_b04.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode05/*.fastq.gz > R37_b05.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode06/*.fastq.gz > R37_b06.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode07/*.fastq.gz > R37_b07.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode08/*.fastq.gz > R37_b08.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode09/*.fastq.gz > R37_b09.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode10/*.fastq.gz > R37_b10.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode11/*.fastq.gz > R37_b11.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode12/*.fastq.gz > R37_b12.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode13/*.fastq.gz > R37_b13.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode14/*.fastq.gz > R37_b14.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode15/*.fastq.gz > R37_b15.fq.gz
cat /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/sup_ON001-DNA-R00037_fast5_pass/pass/barcode16/*.fastq.gz > R37_b16.fq.gz
```

### Perform adapter trimming

```bash
cd ..
# Make sample name list
vi sample_list.txt

"
R08_b01
R08_b02
R08_b03
R08_b04
R08_b05
R08_b06
R37_b01
R37_b02
R37_b03
R37_b04
R37_b05
R37_b06
R37_b07
R37_b08
R37_b09
R37_b10
R37_b11
R37_b12
R37_b13
R37_b14
R37_b15
R37_b16
"

mkdir trimmed
while read f;
    do 
    echo "Trimming sample $f"
    porechop -i ./fastq/$f.fq.gz -o ./trimmed/$f.trim.fq.gz --threads 48 --format fastq.gz > ./trimmed/$f.trim.info.txt
    done < sample_list.txt
```

### Filter read length and map reads to 10-repeat hg38 FLG reference sequence

```bash
while read f;
    do
    echo "Filtering sample $f"
    filtlong --min_length 12000 ./trimmed/$f.trim.fq.gz | gzip > ./trimmed/$f.filt.fastq.gz;
    done < sample_list.txt

mkdir bam_hg38
ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

while read f;
    do
    echo "Mapping sample $f"
    minimap2 -t 12 -ax map-ont --secondary=no $ref ./trimmed/$f.filt.fastq.gz | samtools view -b -o ./bam_hg38/$f.mm-hg38.bam
    done < sample_list.txt
```

### Filter reads by FLAG, QUAL, read length, downsampling to 100X and capturing "full-length" reads (fl_capture.py)

```bash
# Created a bash script containing these commands (downsampling.sh)
"
#!/bin/bash

QUAL=60
LENGTH=12000
SAMFLAG=2308

for file in /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/analysis_1223/bam_hg38/*.mm-hg38.bam; do
filename="${file##*/}"
echo "Processing $file"

# Filter for primary alignments
echo "Filtering for primary alignments with quality score "$QUAL""
samtools view -F $SAMFLAG -q $QUAL -b "$file" > "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.bam"

# Sort BAM and create BED
samtools sort -@ 4 "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.bam" -o "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam"
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam"
bedtools bamtobed -i "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam" > "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bed"

# Capture full-length reads and downsample to 100X (creates a ".fl100sr.bam" file)
python fl_capture.py "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bed" "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam" chr1:152302275-152314808 0.98 100

# Index BAM
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.fl100sr.bam"

# Remove intermediate files
echo "removing intermediate files..."

rm "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.bam"
echo "removed "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.bam""
rm "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam"
echo "removed "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam""
rm "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam.bai"
echo "removed "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bam.bai""
rm "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bed"
echo "removed "${filename%%.*}_F${SAMFLAG}_q${QUAL}_${LENGTH}_headered.sorted.bed""

done
"
chmod u+x downsampling.sh
cd bam_hg38
# Run script
../downsampling.sh
cd ..
```

### Detect repeat copies in reads and separate them (detect_copy.py)

```bash
# Run script
while read f;
    do
        python detect_copies.py ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.bam
        samtools index ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.filt.bam
    done < sample_list.txt
# Create repeat count summary file
truncate -s 0 repeat_count_summary.tsv
while read f;
    do 
        G=`tail -n +2 ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.summmary.tsv | awk -F '\t' '{if ($3 == "yes") print $1}'  | wc -l`
        if [[ "$G" -eq 1 ]]; then
            J=`grep 'yes' ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.summmary.tsv | cut -f 1`
            H=`echo $J | sed 's/-rpt.*//g'`
            echo -e ${f}"\t"${H}"_"${H} >> repeat_count_summary.tsv
        else
            myarr=()
            for i in `grep 'yes' ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.summmary.tsv | cut -f 1`;
                do 
                    H=`echo $i | sed 's/-rpt.*//g'`
                    myarr[${#myarr[@]}]=$H
                done
            
            echo -e ${f}"\t"$((myarr[0]<myarr[1] ? myarr[0] : myarr[1]))"_"$((myarr[0]<myarr[1] ? myarr[1] : myarr[0])) >> repeat_count_summary.tsv;
        fi
    done < sample_list.txt
```

### Map 11/12 repeat reads to the synthetic 11-repeat reference (11-RPT-8.1), whereas 10 repeat reads to hg38

```bash
mkdir bam_alt_ref
while read f;
    do
        if [ -f ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.10rep.fastq ]; then
            echo "Mapping sample $f"
            ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
            minimap2 -t 12 -ax map-ont --secondary=no $ref ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.10rep.fastq | samtools sort - -o ./bam_alt_ref/$f.10rep.hg38.bam
            samtools index ./bam_alt_ref/$f.10rep.hg38.bam
        fi
        if [ -f ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.11rep.fastq ]; then
            echo "Mapping sample $f"
            ref="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_11-RPT-8.1p.fa"
            minimap2 -t 12 -ax map-ont --secondary=no $ref ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.11rep.fastq | samtools sort - -o ./bam_alt_ref/$f.11rep.11ref.bam
            samtools index ./bam_alt_ref/$f.11rep.11ref.bam
        fi
        if [ -f ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.12rep.fastq ]; then
            echo "Mapping sample $f"
            ref="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_11-RPT-8.1p.fa"
            minimap2 -t 12 -ax map-ont --secondary=no $ref ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.12rep.fastq | samtools sort - -o ./bam_alt_ref/$f.12rep.11ref.bam
            samtools index ./bam_alt_ref/$f.12rep.11ref.bam
        fi;
    done < sample_list.txt
```

### Call variants using PMDV on haplotypes separately

```bash
mkdir pmdv
while read f;
    do
        SAMPLE=$f
        echo "Processing PMDV for ${SAMPLE}"
        BASE="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/analysis_1223"
        INPUT_DIR=${BASE}/"bam_alt_ref"
        THREADS="12"
        OUTPUT_DIR=${BASE}"/pmdv/"
        # For 10rep
        if [ -f ./bam_alt_ref/${f}.10rep.hg38.bam ]; then
            OUT_PREFIX=${SAMPLE}"_10rep"
            REF_DIR="/public-data/references/GRCh38"
            REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
            BAM=${f}".10rep.hg38.bam"
            # Run Pepper-margin-deepvariant
            docker run \
            --user "$(id -u):$(id -g)" \
            -v "${INPUT_DIR}":"${INPUT_DIR}" \
            -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
            -v "${REF_DIR}":"${REF_DIR}" \
            kishwars/pepper_deepvariant:r0.8 \
            run_pepper_margin_deepvariant call_variant \
            -b "${INPUT_DIR}/${BAM}" \
            -f "${REF_DIR}/${REF}" \
            -o "${OUTPUT_DIR}" \
            -t "${THREADS}" \
            -p "${OUT_PREFIX}" \
            --ont_r9_guppy5_sup
            if [ ! -f ./pmdv/${f}"_10rep.vcf.gz" ]; then
                zcat ./pmdv/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz | grep -v "refCall" | grep -v '0/1' > ./pmdv/${SAMPLE}"_10rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_10rep.vcf"
            fi
            G=`tail -n +2 ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.summmary.tsv | awk -F '\t' '{if ($3 == "yes") print $1}'  | wc -l`
            if [[ "$G" -eq 1 ]]; then 
                zcat ./pmdv/${SAMPLE}"_10rep.vcf.gz" | grep -v 'refCall' > ./pmdv/${SAMPLE}"_10rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_10rep.vcf"
            else
                zcat ./pmdv/${SAMPLE}"_10rep.vcf.gz" | grep -v "refCall" > ./pmdv/${SAMPLE}"_10rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_10rep.vcf"
            fi
        fi
        # For 11rep
        if [ -f ./bam_alt_ref/${f}.11rep.11ref.bam ]; then
            OUT_PREFIX=${SAMPLE}"_11rep"
            REF_DIR="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/"
            REF="FLG_11-RPT-8.1p.fa"
            BAM=${f}".11rep.11ref.bam"
            # Run Pepper-margin-deepvariant
            docker run \
            --user "$(id -u):$(id -g)" \
            -v "${INPUT_DIR}":"${INPUT_DIR}" \
            -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
            -v "${REF_DIR}":"${REF_DIR}" \
            kishwars/pepper_deepvariant:r0.8 \
            run_pepper_margin_deepvariant call_variant \
            -b "${INPUT_DIR}/${BAM}" \
            -f "${REF_DIR}/${REF}" \
            -o "${OUTPUT_DIR}" \
            -t "${THREADS}" \
            -p "${OUT_PREFIX}" \
            --ont_r9_guppy5_sup
            if [ ! -f ./pmdv/${f}"_11rep.vcf.gz" ]; then
                zcat ./pmdv/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz | grep -v "refCall" | grep -v '0/1' > ./pmdv/${SAMPLE}"_11rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_11rep.vcf"
            fi
            G=`tail -n +2 ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.summmary.tsv | awk -F '\t' '{if ($3 == "yes") print $1}'  | wc -l`
            if [[ "$G" -eq 1 ]]; then 
                zcat ./pmdv/${SAMPLE}"_11rep.vcf.gz" | grep -v 'refCall' > ./pmdv/${SAMPLE}"_11rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_11rep.vcf"
            else
                zcat ./pmdv/${SAMPLE}"_11rep.vcf.gz" | grep -v "refCall" > ./pmdv/${SAMPLE}"_11rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_11rep.vcf"
            fi
        fi
        # For 12rep
        if [ -f ./bam_alt_ref/${f}.12rep.11ref.bam ]; then
            OUT_PREFIX=${SAMPLE}"_12rep"
            REF_DIR="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/"
            REF="FLG_11-RPT-8.1p.fa"
            BAM=${f}".12rep.11ref.bam"
            # Run Pepper-margin-deepvariant
            docker run \
            --user "$(id -u):$(id -g)" \
            -v "${INPUT_DIR}":"${INPUT_DIR}" \
            -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
            -v "${REF_DIR}":"${REF_DIR}" \
            kishwars/pepper_deepvariant:r0.8 \
            run_pepper_margin_deepvariant call_variant \
            -b "${INPUT_DIR}/${BAM}" \
            -f "${REF_DIR}/${REF}" \
            -o "${OUTPUT_DIR}" \
            -t "${THREADS}" \
            -p "${OUT_PREFIX}" \
            --ont_r9_guppy5_sup
            if [ ! -f ./pmdv/${f}"_12rep.vcf.gz" ]; then
                zcat ./pmdv/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz | grep -v "refCall" | grep -v '0/1' > ./pmdv/${SAMPLE}"_12rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_12rep.vcf"
            fi
            G=`tail -n +2 ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.summmary.tsv | awk -F '\t' '{if ($3 == "yes") print $1}'  | wc -l`
            if [[ "$G" -eq 1 ]]; then 
                zcat ./pmdv/${SAMPLE}"_12rep.vcf.gz" | grep -v 'refCall' > ./pmdv/${SAMPLE}"_12rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_12rep.vcf"
            else
                zcat ./pmdv/${SAMPLE}"_12rep.vcf.gz" | grep -v "refCall" > ./pmdv/${SAMPLE}"_12rep.vcf"
                bgzip -f ./pmdv/${SAMPLE}"_12rep.vcf"
            fi
        fi
    done < sample_list.txt
```

### Liftover coordinates in PMDV 11 ref VCF to match hg38 (Awk)

```bash
mkdir liftover
while read f;
    do
        # 10rep
        if [ -f ./pmdv/${f}"_10rep.vcf.gz" ]; then
            cp ./pmdv/${f}"_10rep.vcf.gz" ./liftover/${f}_10rep.vcf.gz
            tabix ./liftover/${f}_10rep.vcf.gz
        fi
        # 11rep
        if [ -f ./pmdv/${f}"_11rep.vcf.gz" ]; then
            bcftools view -h ./pmdv/${f}"_11rep.vcf.gz" > temp_header.txt
            bcftools view -H ./pmdv/${f}"_11rep.vcf.gz" | awk -F'\t' '{ if ($2 < 6696) $2=$2+152300000; else if ($2 >= 6696) $2=$2+152299028} 1' OFS='\t' > temp_vcf.txt
            cat temp_header.txt temp_vcf.txt | bcftools sort - | perl -pe 's/FLG_11-RPT-8.1p/chr1/g' > temp_vcf2.txt
            bcftools norm -d both temp_vcf2.txt -o ./liftover/${f}_11rep.vcf
            bgzip -f ./liftover/${f}_11rep.vcf
            tabix ./liftover/${f}_11rep.vcf.gz
            rm temp_header.txt
            rm temp_vcf.txt
            rm temp_vcf2.txt
        fi
        # 12rep
        if [ -f ./pmdv/${f}"_12rep.vcf.gz" ]; then
            bcftools view -h ./pmdv/${f}"_12rep.vcf.gz" > temp_header.txt
            bcftools view -H ./pmdv/${f}"_12rep.vcf.gz" | awk -F'\t' '{ if ($2 < 6696) $2=$2+152300000; else if ($2 >= 6696) $2=$2+152299028} 1' OFS='\t' > temp_vcf.txt
            cat temp_header.txt temp_vcf.txt | bcftools sort - | perl -pe 's/FLG_11-RPT-8.1p/chr1/g' > temp_vcf2.txt
            bcftools norm -d both temp_vcf2.txt -o ./liftover/${f}_12rep.vcf
            bgzip -f ./liftover/${f}_12rep.vcf
            tabix ./liftover/${f}_12rep.vcf.gz
            rm temp_header.txt
            rm temp_vcf.txt
            rm temp_vcf2.txt
        fi
    done < sample_list.txt
```

### Merge haplotypes using custom script (haplo_merge.py) and haplotag BAM (Whatshap)

```bash
mkdir merged
while read f;
    do
        C=`ls ./liftover/${f}*rep.vcf.gz | wc -l`
        if [[ "$C" -eq 2 ]]; then  # If two haplotypes found
            myarr=()
            for i in ./liftover/${f}*rep.vcf.gz;
                do 
                    myarr[${#myarr[@]}]=$i
                done
            python haplo_merge.py ${myarr[0]} ${myarr[1]} | bcftools sort - > ./merged/${f}_merged.vcf
            bgzip -f ./merged/${f}_merged.vcf
        elif [[ "$C" -eq 1 ]]; then
            for i in ./liftover/${f}*rep.vcf.gz;
                do
                    cp $i ./merged/${f}_merged.vcf.gz
                done
        fi
        tabix ./merged/${f}_merged.vcf.gz
    done < sample_list.txt
```

### Phase variants with Whatshap

```bash
ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
mkdir whatshap
while read f;
    do
        echo "Running WhatsHap on sample $f"
        whatshap phase -o ./whatshap/${f}_phased.vcf --ignore-read-groups --indels --reference=$ref \
        ./merged/${f}_merged.vcf.gz \
        ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.filt.bam
        bgzip --keep ./whatshap/${f}_phased.vcf
        tabix ./whatshap/${f}_phased.vcf.gz
        whatshap haplotag -o ./whatshap/${f}.haplotagged.bam --ignore-read-groups --reference=$ref \
        ./whatshap/${f}_phased.vcf.gz \
        ./bam_hg38/${f}_F2308_q60_12000_headered.sorted.fl100sr.filt.bam
        samtools index ./whatshap/${f}.haplotagged.bam;
    done < sample_list.txt
```

### Perform known variant check

```bash
# Parse variant of interest file to individual sample files
mkdir variants_of_interest && cd variants_of_interest
tail -n +3 variants.csv  | cut -d',' -f 2,5,11 | perl -pe 's/,,\n//g' | perl -pe 's/R000/R/g' | perl -pe 's/barcode/b/g' | perl -pe 's/,NIL//g' > variants2.csv
awk -F, '{print $2"\n"$3 > $1".txt"}' variants2.csv
cd ..
# Remove empty line in R08_b01

# Check presence of variant coordinates in VCFs
while read f;
    do 
        echo "----------- Checking sample $f -----------"
        echo "------ PMDV $f ------"
        awk -F'\t' 'NR==FNR{c[$1]++;next};$2 in c' ./variants_of_interest/${f}.txt ./whatshap/${f}_phased.vcf
        echo "-"
        echo "-";
    done < sample_list.txt
```
