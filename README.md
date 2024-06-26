# Repository for "Nanopore Sequencing Enables Allelic Phasing of FLG Loss-of-Function Variants, Intragenic Copy Number Variation, and Methylation Status in Atopic Dermatitis and Ichthyosis Vulgaris"

## Publication
Wong C, Tham CY, Yang L, Benton MC, Narang V, Denil S, Duan K, Yew YW, Lee B, Florez de Sessions P, Common JEA. Nanopore Sequencing Enables Allelic Phasing of FLG Loss-of-Function Variants, Intragenic Copy Number Variation, and Methylation Status in Atopic Dermatitis and Ichthyosis Vulgaris. J Invest Dermatol. 2024 Feb 8:S0022-202X(24)00097-6. doi: [10.1016/j.jid.2024.01.020](https://doi.org/10.1016/j.jid.2024.01.020). PMID: 38336337.

## Contents
| File | Description |
| --- | --- |
| [amplicon_analysis_pipeline.md](./scripts/amplicon_analysis_pipeline.md) | Code for FLG amplicon analysis |
| [adaptive-sampling_analysis_pipeline.md](./scripts/adaptive-sampling_analysis_pipeline.md) | Code for FLG adaptive sampling analysis |
| [detect_copies.py](./scripts/detect_copies.py) | Python script for detecting FLG RPT copies in amplicon reads |
| [detect_copies_AS.py](./scripts/detect_copies_AS.py) | Python script for detecting FLG RPT copies in adaptive sampling reads|
| [fl_capture.py](./scripts/fl_capture.py) | Python script to capture close to full-length FLG amplicon reads |
| [haplo_merge.py](./scripts/haplo_merge.py) | Python script to merge VCF files from two haplotypes |
| [FLG_11-RPT-8.1p.fa](./scripts/FLG_11-RPT-8.1p.fa) | Synthetic FLG reference 11-RPT-8.1' FASTA |
| [FLG_12-RPT-8.1p-10.1p.fa](./scripts/FLG_12-RPT-8.1p-10.1p.fa) | Synthetic FLG reference 12-RPT-8.1'-10.1' FASTA |
