# Filaggrin adaptive sampling variant calling, phasing, and modbasecalling bioinformatics pipeline

## Pipeline to generate phased BAMs for methylation visualization

### Sample ID legend

```python
sample_legend = {
    '3A': 'S11_11', 
    '6A': 'S12_12' 
}
```

### Modbasecalling

```bash
# S11_11
GUPPYPATH="/public-data/software/guppy/6.4.2/ont-guppy/bin/"
REF="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
INPUTDIR="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11"
OUTPUTDIR="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall"
CONFIG="dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_sup.cfg"
"${GUPPYPATH}"/guppy_basecaller \
    -i "${INPUTDIR}" \
    -r \
    -s "${OUTPUTDIR}" \
    -x 'cuda:0,1' \
    -c "${CONFIG}" \
    -a "${REF}" \
    --chunks_per_runner 412 \
    --bam_out

# Merge BAM files and index
samtools cat -@ 24 -o S11_11_pass.bam pass/*.bam
samtools sort -@ 24 S11_11_pass.bam -o S11_11_pass.sort.bam
samtools index S11_11_pass.sort.bam
rm S11_11_pass.bam

# S12_12
GUPPYPATH="/public-data/software/guppy/6.4.2/ont-guppy/bin/"
REF="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
INPUTDIR="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12"
OUTPUTDIR="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall"
CONFIG="dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_sup.cfg"
"${GUPPYPATH}"/guppy_basecaller \
    -i "${INPUTDIR}" \
    -r \
    -s "${OUTPUTDIR}" \
    -x 'cuda:0,1' \
    -c "${CONFIG}" \
    -a "${REF}" \
    --chunks_per_runner 412 \
    --bam_out

# Merge BAM files and index
samtools cat -@ 24 -o S12_12_pass.bam pass/*.bam
samtools sort -@ 24 S12_12_pass.bam -o S12_12_pass.sort.bam
samtools index S12_12_pass.sort.bam
rm S12_12_pass.bam
```

### Filter reads by FLAG, QUAL, and mapping region (+-20kb of FLG gene), and detect FLG RPT copies

```bash
# S11_11
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall
mkdir analysis && cd analysis

file="../S11_11_pass.sort.bam"
filename="${file##*/}"
QUAL=60
SAMFLAG=2308
samtools view -@ 24 -F $SAMFLAG -q $QUAL -b "$file" > "${filename%%.*}_F${SAMFLAG}_q${QUAL}.bam"
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}.bam"
samtools view -@ 24 -F $SAMFLAG -q $QUAL -b "$file" chr1:152302165-152325239 > "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"  # range used was +-20kb of FLG gene
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"

python detect_copies_AS.py "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"
samtools index ./bam_hg38/${f}_F2308_q60_12000_headered.fl100sr.filt.bam

# S12_12
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall
mkdir analysis && cd analysis

file="../S12_12_pass.sort.bam"
filename="${file##*/}"
QUAL=60
SAMFLAG=2308
samtools view -@ 24 -F $SAMFLAG -q $QUAL -b "$file" > "${filename%%.*}_F${SAMFLAG}_q${QUAL}.bam"
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}.bam"
samtools view -@ 24 -F $SAMFLAG -q $QUAL -b "$file" chr1:152292165-152335239 > "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"  # range used was +-20kb of FLG gene
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"

python detect_copies_AS.py "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"
samtools index ./bam_hg38/${f}_F2308_q60_12000_headered.fl100sr.filt.bam
```

### Call variants using PMDV

```bash
# S11_11
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis
samtools fastq S11_11_pass_F2308_q60_FLG.bam > S11_11_pass_F2308_q60_FLG.fq
REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_11-RPT-8.1p.fa'
minimap2 -ax map-ont $REF S11_11_pass_F2308_q60_FLG.fq | samtools sort -o S11_11_pass_F2308_q60_FLG-11ref.bam
samtools index S11_11_pass_F2308_q60_FLG-11ref.bam
mkdir pmdv
SAMPLE="S11_11"
echo "Processing PMDV for ${SAMPLE}"
BASE="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/${SAMPLE}/modbasecall/analysis"
INPUT_DIR="${BASE}"
THREADS="12"
OUTPUT_DIR=${BASE}"/pmdv/"
OUT_PREFIX=${SAMPLE}_pmdv
REF_DIR="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5"
REF="FLG_11-RPT-8.1p.fa"
BAM=${SAMPLE}"_pass_F2308_q60_FLG-11ref.bam"
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

# S12_12
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis
samtools fastq S12_12_pass_F2308_q60_FLG.bam > S12_12_pass_F2308_q60_FLG.fq
REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_11-RPT-8.1p.fa'
minimap2 -ax map-ont $REF S12_12_pass_F2308_q60_FLG.fq | samtools sort -o S12_12_pass_F2308_q60_FLG-11ref.bam
samtools index S12_12_pass_F2308_q60_FLG-11ref.bam
mkdir pmdv
SAMPLE="S12_12"
echo "Processing PMDV for ${SAMPLE}"
BASE="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/${SAMPLE}/modbasecall/analysis"
INPUT_DIR="${BASE}"
THREADS="12"
OUTPUT_DIR=${BASE}"/pmdv/"
OUT_PREFIX=${SAMPLE}_pmdv
REF_DIR="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5"
REF="FLG_11-RPT-8.1p.fa"
BAM=${SAMPLE}"_pass_F2308_q60_FLG-11ref.bam"
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
```

### Liftover coordinates in PMDV 11 ref VCF to match hg38 (Awk)

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis
f='S11_11'
mkdir liftover
bcftools view -h ./pmdv/${f}"_pmdv.vcf.gz" > temp_header.txt
bcftools view -H ./pmdv/${f}"_pmdv.vcf.gz" | awk -F'\t' '{ if ($2 < 6696) $2=$2+152300000; else if ($2 >= 6696) $2=$2+152299028} 1' OFS='\t'  > temp_vcf.txt
cat temp_header.txt temp_vcf.txt | bcftools sort - | perl -pe 's/FLG_11-RPT-8.1p/chr1/g' > temp_vcf2.txt
bcftools norm -d both temp_vcf2.txt -o ./liftover/${f}_pmdv_lift.vcf
bgzip -f ./liftover/${f}_pmdv_lift.vcf
tabix ./liftover/${f}_pmdv_lift.vcf.gz
rm temp_header.txt
rm temp_vcf.txt
rm temp_vcf2.txt

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis
f='S12_12'
mkdir liftover
bcftools view -h ./pmdv/${f}"_pmdv.vcf.gz" > temp_header.txt
bcftools view -H ./pmdv/${f}"_pmdv.vcf.gz" | awk -F'\t' '{ if ($2 < 6696) $2=$2+152300000; else if ($2 >= 6696) $2=$2+152299028} 1' OFS='\t'  > temp_vcf.txt
cat temp_header.txt temp_vcf.txt | bcftools sort - | perl -pe 's/FLG_11-RPT-8.1p/chr1/g' > temp_vcf2.txt
bcftools norm -d both temp_vcf2.txt -o ./liftover/${f}_pmdv_lift.vcf
bgzip -f ./liftover/${f}_pmdv_lift.vcf
tabix ./liftover/${f}_pmdv_lift.vcf.gz
rm temp_header.txt
rm temp_vcf.txt
rm temp_vcf2.txt
```

### Phase variants with Whatshap

```bash
ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis
mkdir whatshap
f="S11_11"
echo "Running WhatsHap on sample $f"
whatshap phase -o ./whatshap/${f}_phased.vcf --ignore-read-groups --indels --reference=$ref \
./liftover/${f}_pmdv_lift.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
bgzip --keep ./whatshap/${f}_phased.vcf
tabix ./whatshap/${f}_phased.vcf.gz
whatshap haplotag -o ./whatshap/${f}.FLG.haplotagged.bam --ignore-read-groups \
--reference=$ref \
./whatshap/${f}_phased.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
samtools index ./whatshap/${f}.FLG.haplotagged.bam;
# For full genome
whatshap haplotag -o ./whatshap/${f}.haplotagged.bam --ignore-read-groups \
--reference=$ref \
--skip-missing-contigs \
./whatshap/${f}_phased.vcf.gz \
./${f}_pass_F2308_q60.bam
samtools index ./whatshap/${f}.haplotagged.bam;

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis
mkdir whatshap
f="S12_12"
echo "Running WhatsHap on sample $f"
whatshap phase -o ./whatshap/${f}_phased.vcf --ignore-read-groups --indels --reference=$ref \
./liftover/${f}_pmdv_lift.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
bgzip --keep ./whatshap/${f}_phased.vcf
tabix ./whatshap/${f}_phased.vcf.gz
whatshap haplotag -o ./whatshap/${f}.FLG.haplotagged.bam --ignore-read-groups \
--reference=$ref \
./whatshap/${f}_phased.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
samtools index ./whatshap/${f}.FLG.haplotagged.bam;
# For full genome
whatshap haplotag -o ./whatshap/${f}.haplotagged.bam --ignore-read-groups \
--reference=$ref \
--skip-missing-contigs \
./whatshap/${f}_phased.vcf.gz \
./${f}_pass_F2308_q60.bam
samtools index ./whatshap/${f}.haplotagged.bam;
```

## Pipeline to generate phased BAMs with extended promoter region for methylation visualization

### Filter reads by FLAG, QUAL, and mapping region (Extended)

```bash
# S11_11
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall
mkdir analysis_extended && cd analysis_extended

file="../S11_11_pass.sort.bam"
filename="${file##*/}"
QUAL=60
SAMFLAG=2308
samtools view -@ 24 -F $SAMFLAG -q $QUAL -b "$file" chr1:152090000-152530000 > "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"

# S12_12
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall
mkdir analysis_extended && cd analysis_extended

file="../S12_12_pass.sort.bam"
filename="${file##*/}"
QUAL=60
SAMFLAG=2308
samtools view -@ 24 -F $SAMFLAG -q $QUAL -b "$file" chr1:152090000-152530000 > "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"
samtools index "${filename%%.*}_F${SAMFLAG}_q${QUAL}_FLG.bam"
```

### Extend FLG_11-RPT-8.1p reference fasta

```bash
cd /apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5

# Original start
chr1:152300001

# Original end
chr1:152317000

# Sequences to add
chr1:152090000-152300000
chr1:152317000-152530000

seqtk subseq /public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna FLG_extended_coord.bed > FLG_extended_seq.fa

grep 'chr1:152090001'  FLG_extended_seq.fa -A1 | tail -n +2 | tr -d '\n' > s1.tmp
grep 'chr1:152317001'  FLG_extended_seq.fa -A1 | tail -n +2 > s3.tmp
tail -n +2 FLG_11-RPT-8.1p.fa | tr -d '\n' > s2.tmp
cat s1.tmp s2.tmp s3.tmp > s4.tmp 
echo ">FLG_11-RPT-8.1p_ext" | cat - s4.tmp > FLG_11-RPT-8.1p_ext.fa
samtools faidx FLG_11-RPT-8.1p_ext.fa
rm *.tmp
```

### Call variants using PMDV against FLG_11-RPT-8.1p_ext reference

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis_extended
samtools fastq S11_11_pass_F2308_q60_FLG.bam > S11_11_pass_F2308_q60_FLG.fq
REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_11-RPT-8.1p_ext.fa'
minimap2 -ax map-ont $REF S11_11_pass_F2308_q60_FLG.fq | samtools sort -o S11_11_pass_F2308_q60_FLG-11ref.bam
samtools index S11_11_pass_F2308_q60_FLG-11ref.bam
mkdir pmdv
SAMPLE="S11_11"
echo "Processing PMDV for ${SAMPLE}"
BASE="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/${SAMPLE}/modbasecall/analysis_extended"
INPUT_DIR="${BASE}"
THREADS="12"
OUTPUT_DIR=${BASE}"/pmdv/"
OUT_PREFIX=${SAMPLE}_pmdv
REF_DIR="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5"
REF="FLG_11-RPT-8.1p_ext.fa"
BAM=${SAMPLE}"_pass_F2308_q60_FLG-11ref.bam"
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

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis_extended
samtools fastq S12_12_pass_F2308_q60_FLG.bam > S12_12_pass_F2308_q60_FLG.fq
REF='/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5/FLG_11-RPT-8.1p_ext.fa'
minimap2 -ax map-ont $REF S12_12_pass_F2308_q60_FLG.fq | samtools sort -o S12_12_pass_F2308_q60_FLG-11ref.bam
samtools index S12_12_pass_F2308_q60_FLG-11ref.bam
mkdir pmdv
SAMPLE="S12_12"
echo "Processing PMDV for ${SAMPLE}"
BASE="/apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/${SAMPLE}/modbasecall/analysis_extended"
INPUT_DIR="${BASE}"
THREADS="12"
OUTPUT_DIR=${BASE}"/pmdv/"
OUT_PREFIX=${SAMPLE}_pmdv
REF_DIR="/apps-data/202202_filaggrin_AS/Nanopore_FLG_amplicon_seq_fast5"
REF="FLG_11-RPT-8.1p_ext.fa"
BAM=${SAMPLE}"_pass_F2308_q60_FLG-11ref.bam"
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
```

### Liftover coordinates in PMDV 11 ref VCF to match hg38 (Awk) (Extended)

```bash
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis_extended
f='S11_11'
mkdir liftover
bcftools view -h ./pmdv/${f}"_pmdv.vcf.gz" > temp_header.txt
bcftools view -H ./pmdv/${f}"_pmdv.vcf.gz" | awk -F'\t' '{ if ($2 < 216498) $2=$2+152090000; else if ($2 >= 216498) $2=$2+152089028} 1' OFS='\t'  > temp_vcf.txt
cat temp_header.txt temp_vcf.txt | bcftools sort - | perl -pe 's/FLG_11-RPT-8.1p_ext/chr1/g' > ./temp_vcf2.txt
bcftools norm -d both temp_vcf2.txt -o ./liftover/${f}_11rep.vcf
bgzip -f ./liftover/${f}_11rep.vcf
tabix ./liftover/${f}_11rep.vcf.gz
rm temp_header.txt
rm temp_vcf.txt
rm temp_vcf2.txt

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis_extended
f='S12_12'
mkdir liftover
bcftools view -h ./pmdv/${f}"_pmdv.vcf.gz" > temp_header.txt
bcftools view -H ./pmdv/${f}"_pmdv.vcf.gz" | awk -F'\t' '{ if ($2 < 216498) $2=$2+152090000; else if ($2 >= 216498) $2=$2+152089028} 1' OFS='\t'  > temp_vcf.txt
cat temp_header.txt temp_vcf.txt | bcftools sort - | perl -pe 's/FLG_11-RPT-8.1p_ext/chr1/g' > ./temp_vcf2.txt
bcftools norm -d both temp_vcf2.txt -o ./liftover/${f}_11rep.vcf
bgzip -f ./liftover/${f}_11rep.vcf
tabix ./liftover/${f}_11rep.vcf.gz
rm temp_header.txt
rm temp_vcf.txt
rm temp_vcf2.txt
```

### Phase variants with Whatshap (Extended)

```bash
ref="/public-data/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S11_11/modbasecall/analysis_extended
mkdir whatshap
f="S11_11"
echo "Running WhatsHap on sample $f"
whatshap phase -o ./whatshap/${f}_phased.vcf --ignore-read-groups --indels --reference=$ref \
./liftover/${f}_11rep.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
bgzip --keep ./whatshap/${f}_phased.vcf
tabix ./whatshap/${f}_phased.vcf.gz
whatshap haplotag -o ./whatshap/${f}.FLG.haplotagged_ext.bam --ignore-read-groups \
--reference=$ref \
./whatshap/${f}_phased.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
samtools index ./whatshap/${f}.FLG.haplotagged_ext.bam;

cd /apps-data/202202_filaggrin_AS/FLG_AS/EDC_AS_11May/S12_12/modbasecall/analysis_extended
mkdir whatshap
f="S12_12"
echo "Running WhatsHap on sample $f"
whatshap phase -o ./whatshap/${f}_phased.vcf --ignore-read-groups --indels --reference=$ref \
./liftover/${f}_11rep.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
bgzip --keep ./whatshap/${f}_phased.vcf
tabix ./whatshap/${f}_phased.vcf.gz
whatshap haplotag -o ./whatshap/${f}.FLG.haplotagged_ext.bam --ignore-read-groups \
--reference=$ref \
./whatshap/${f}_phased.vcf.gz \
./${f}_pass_F2308_q60_FLG.bam
samtools index ./whatshap/${f}.FLG.haplotagged_ext.bam;
```

## FLG adaptive sampling methylation analysis pipeline

### Get FLG gene feature coordinates from GTF file

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

### Run modkit for the various FLG gene features

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

### Calculate methylation percentage for gene features

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

### Map uBAM to 11/12 rep reference for methylation analysis on each RPT regions

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

### Get RPT coordinates from Sandilands et al., 2007 (10.1038/ng2020)

```bash
"""
# Refer to FLG_RPT_coord_sandilands.xlsx (Not provided) (Inferred from Sandilands et al., 2007, sup fig 1)

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

## Run modkit for each RPT region

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

## Calculate methylation percentage for each RPT region

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

# Checking for methylation percentage
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

# Checking for methylation percentage
for i in *.bed; do
    echo ${i}
    awk '{if (($4=="m") && ($6=="+")) print $0}' ${i} | cut -f 10 | cut -f 1,4 -d ' ' | awk '{s+=$1}{j+=$2}END{print (1-(j/s))*100}'
    echo ""
    awk '{if (($4=="m") && ($6=="-")) print $0}' ${i} | cut -f 10 | cut -f 1,4 -d ' ' | awk '{s+=$1}{j+=$2}END{print (1-(j/s))*100}'
    echo ""
done
```
