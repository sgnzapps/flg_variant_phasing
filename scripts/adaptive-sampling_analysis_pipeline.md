# Filaggrin adaptive sampling variant calling, phasing, and modbasecalling bioinformatics pipeline

## Pipeline

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
