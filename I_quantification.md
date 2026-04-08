---
output:
  html_document: default
  pdf_document: default
---

# Natural variation in transplacental transfer efficiency exposes distinct transcriptional network architectures of PFAS effects on birth weight and gestational age

## I. Short-read quantification workflow

**Overview**: This section presents a comprehensive analysis of short-read RNA-seq data from the GUSTO birth cohort and independent cohort of experimental patient-derived placental explants to explore the transcriptional network architectures underlying PFAS effects on perinatal outcomes.

**Note**: The placenta reference transcriptome (lr-assembly) can be downloaded from our GitHub repository, [lr-placenta-transcriptomics/Assembly/ESPRESSO_corrected_SC_filtered.gtf.gz](https://github.com/sbresnahan/lr-placenta-transcriptomics/AssemblyESPRESSO_corrected_SC_filtered.gtf.gz)

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Primary data directories
DIR_RAW=/path/to/raw_reads                    # Raw sequencing data
DIR_TRIM=/path/to/trimmed_reads               # Quality-trimmed reads
DIR_COUNTS_ASSEMBLY=/path/to/assembly_counts  # Salmon quantification (assembly)
DIR_COUNTS_GENCODE=/path/to/gencode_counts    # Salmon quantification (GENCODE)

# Processing and temporary directories
DIR_INDEX=/path/to/preprocessing_indices_directory    # Reference genome indices
DIR_REPORTS=/path/to/reports_directory                # QC and processing reports
DIR_METRICS=/path/to/metrics_directory                # Alignment and trimming metrics

# Reference genome and annotation files
GENOME_FASTA=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna  # hg38 genome
REFERENCE_GTF=/path/to/gencode_v45_annotation.gtf                      # GENCODE v45
```

#### Generate Salmon mapping index for lr-assembly

**Purpose**: Create decoy-aware Salmon index for the filtered placental transcriptome assembly. Decoy sequences from the genome help distinguish between transcriptomic and genomic origins of reads, improving quantification accuracy and reducing false positive assignments.

```bash
################################################################################
# Directory Setup and Configuration
################################################################################

# Reference genome directory
DIR_GENOME=/path/to/genome/directory

# lr-assembly directory
DIR_TXOME=/path/to/lr_assembly/directory
gunzip ${DIR_TXOME}/ESPRESSO_corrected_SC_filtered.gtf.gz

# Output directory for Salmon index
DIR_OUT=${DIR_TXOME}/salmon_index
mkdir -p ${DIR_OUT}

################################################################################
# Generate Decoy-Aware Transcriptome (Gentrome)
################################################################################

# Purpose: Create gentrome (transcriptome + genomic decoys) to improve mapping accuracy
# Decoy sequences help identify reads that might map to genomic regions not 
# represented in the transcriptome, reducing false positive transcript quantification

# Requirements:
#   - generateDecoyTranscriptome.sh (SalmonTools)
#   - bedtools
#   - salmon v1.10.2

# Generate gentrome using SalmonTools generateDecoyTranscriptome.sh
# This script extracts genomic sequences that could act as decoys for transcript mapping
sh generateDecoyTranscriptome.sh \
 -j 10 \
 -a ${DIR_TXOME}/ESPRESSO_corrected_SC_filtered.gtf \
 -g ${GENOME_FASTA} \
 -t ${DIR_TXOME}/ESPRESSO_corrected_SC_filtered.fasta \
 -o ${DIR_OUT}

# Parameter explanations:
# -j 10: Use 10 threads for parallel processing
# -b: Path to bedtools binary for genomic interval operations
# -a: Assembly GTF file with transcript coordinates  
# -g: Reference genome FASTA file
# -t: Transcript sequences FASTA file
# -o: Output directory for gentrome files

################################################################################
# Build Salmon Index
################################################################################

# Purpose: Create Salmon k-mer index for fast and accurate transcript quantification
# Uses decoy-aware indexing to account for genomic sequences that could 
# interfere with transcript mapping

# Build Salmon index with decoy awareness
salmon index -t ${DIR_OUT}/gentrome.fa \
  -i ${DIR_OUT}/ESPRESSO_corrected_SC_filtered \
  --decoys ${DIR_OUT}/decoys.txt \
  -k 31
  
################################################################################
# Technical Notes
################################################################################

# Parameter explanations:
# -t: Gentrome FASTA file (transcripts + decoy sequences)
# -i: Output index directory name
# --decoys: Text file listing decoy sequence names
# -k 31: k-mer size (31 is optimal for most RNA-seq applications)

# Decoy-aware indexing benefits:
# 1. Reduces false positive mappings to transcripts
# 2. Accounts for genomic regions not represented in transcriptome
# 3. Improves quantification accuracy, especially for lowly expressed transcripts
# 4. Handles reads from unannotated genomic regions

# k-mer size considerations:
# k=31 provides good balance between:
# - Specificity (longer k-mers are more specific)
# - Sensitivity (shorter k-mers capture more mappings)
# - Index size (longer k-mers create larger indices)

# Expected index characteristics:
# - Size: ~2-4 GB for human transcriptome + decoys
# - Build time: ~30-60 minutes depending on hardware
# - Memory usage: ~8-16 GB during construction
```

#### Generate Salmon mapping index for GENCODEv45

**Purpose**: Create equivalent Salmon index for GENCODE v45 reference to enable direct comparison with assembly quantification. Uses identical parameters and decoy strategy for fair comparison.

```bash
################################################################################
# Directory Setup and Configuration
################################################################################

# GENCODE reference directory
DIR_TXOME_GENCODE=/path/to/GENCODE_v45

# Output directory for Salmon index
DIR_OUT=${DIR_TXOME_GENCODE}/gencode.v45.salmon_index
mkdir -p ${DIR_OUT}

################################################################################
# Generate Decoy-Aware Transcriptome (Gentrome) for GENCODE v45
################################################################################

# Generate gentrome using identical parameters as assembly index
sh generateDecoyTranscriptome.sh \
 -j 10 \
 -a ${REFERENCE_GTF} \
 -g ${GENOME_FASTA} \
 -t ${DIR_TXOME_GENCODE}/gencode.v45.transcripts.fa \
 -o ${DIR_OUT}

################################################################################
# Build Salmon Index for GENCODE v45
################################################################################

# Build Salmon index with identical parameters as assembly
salmon index -t ${DIR_OUT}/gentrome.fa \
  -i ${DIR_OUT}/gencode_v45 \
  --decoys ${DIR_OUT}/decoys.txt \
  -k 31

# Consistent indexing ensures fair comparison between assembly and reference
```

#### Process GUSTO & Explants short reads

**Purpose**: Process short-read RNA-seq data from two independent cohorts for dual quantification against both lr-assembly and reference transcriptomes. This provides validation across multiple datasets and population backgrounds.

**Cohort characteristics**:

- **GUSTO**: Singaporean birth cohort, term placental samples (N=124) with p=8 PFAS measures in maternal blood and cord blood at delivery
- **Explants**: Patient-derived placental explants (N=18) from 2nd trimester and term, treated with PFBS at 0uM, 5uM, and 20uM for 24hrs (n=3 per group)

```bash
################################################################################
# Sample Processing Configuration
################################################################################

# Define sample identifiers for batch processing
LIBIDs=() # SampleIDs as in .fastq file names, e.g. Sample01.fastq.gz

# Processing parameters
THREADS=10  # Parallel processing threads
```

##### Adapter trimming with fastp

**Purpose**: Remove sequencing adapters and low-quality bases to improve mapping accuracy and reduce technical artifacts in quantification.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Quality control and adapter trimming with fastp
# More stringent parameters than typical to ensure high-quality reads
fastp \
  -i ${DIR_RAW}/${LIBID}_1.fq.gz \
  -o ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -I ${DIR_RAW}/${LIBID}_2.fq.gz \
  -O ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -w ${THREADS} \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 36 \
  --html ${DIR_REPORTS}/${LIBID}_fastp.html \
  --json ${DIR_REPORTS}/${LIBID}_fastp.json

# Parameter explanations:
# --cut_tail: Remove low quality bases from 3' end
# --cut_window_size 4: Sliding window size for quality assessment
# --cut_mean_quality 20: Minimum average quality in sliding window
# --length_required 36: Minimum read length after trimming
# --html/--json: Generate detailed QC reports
```

##### Quantification against lr-assembly

**Purpose**: Quantify transcript expression using the filtered placental assembly. Includes inferential replicates for uncertainty assessment and bias correction for improved accuracy.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Salmon quantification against placental assembly

salmon quant \
  -i ${DIR_TXOME}/salmon_index \
  --libType A \
  --validateMappings \
  --seqBias \
  -p ${THREADS} \
  --numBootstraps 100 \
  -1 ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -2 ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -o ${DIR_COUNTS_ASSEMBLY}/${LIBID}

# Parameter explanations:
# --libType A: Automatic library type detection
# --validateMappings: Use more sensitive mapping validation
# --seqBias: Correct for sequence-specific bias
# --numBootstraps 100: Generate inferential replicates for uncertainty
```

##### Quantification against GENCODEv45

**Purpose**: Quantify the same samples against GENCODE v45 reference using identical parameters for direct comparison with assembly results.

For each \$LIBID *n* in \$\{LIBIDs[1...*n*]\}
```bash
# Salmon quantification against GENCODE v45 reference

salmon quant \
  -i ${DIR_TXOME_GENCODE}/gencode.v45.salmon_index \
  --libType A \
  --validateMappings \
  --seqBias \
  -p ${THREADS} \
  --numBootstraps 100 \
  -1 ${DIR_TRIM}/${LIBID}_1.fastq.gz \
  -2 ${DIR_TRIM}/${LIBID}_2.fastq.gz \
  -o ${DIR_COUNTS_GENCODE}/${LIBID}

# Identical parameters ensure fair comparison between annotations
```

---