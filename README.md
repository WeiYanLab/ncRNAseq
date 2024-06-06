# lncRNA and sncRNA Analysis published in journal paper "Epigenetic Transgenerational Inheritance of Toxicant Exposure Specific non-coding RNA in Sperm" - doi:

This repository contains scripts and data for the analysis of long non-coding RNAs (lncRNAs) and small non-coding RNAs (sncRNAs) differential expression using DESeq2.

## Table of Contents

- [Introduction](#introduction)
- [LncRNA Analysis](#lncrna-analysis)
  - [Pipeline Overview](#pipeline-overview-lncrna)
  - [Results](#results-lncrna)
- [sncRNA Analysis](#sncrna-analysis)
  - [Sample Information](#sample-information-sncrna)
  - [Pipeline Overview](#pipeline-overview-sncrna)
  - [Results](#results-sncrna)
- [Scripts and Reference Data](#scripts-and-reference-data)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [License](#license)

## Introduction

This project aims to identify differentially expressed lncRNAs and sncRNAs in various sample groups exposed to different substances (Antrazine, DDT, JetFuel, Vinclozolin) compared to a control group.

## LncRNA Analysis

### Pipeline Overview <a name="pipeline-overview-lncrna"></a>

1. **Alignment**: Hisat2 was used as the aligner.
2. **Counting**: FeatureCounts with parameters `-M` and `-O` for sensitivity.
3. **Reference**: Rattus norvegicus genome (Rnor_6.0) and GTF files from Ensembl [May 2021 archive](https://may2021.archive.ensembl.org/Rattus_norvegicus/Info/Index).
4. **Filtering**: 3249 transcripts filtered from GTF, 29 filtered out using CPAT (coding probability < 0.35).
5. **Differential Expression Analysis**: DESeq2 used with lncRNAs having at least 10 counts. Adjusted p-value < 0.1 and fold change > 2 were thresholds.

## sncRNA Analysis

### Pipeline Overview <a name="pipeline-overview-sncrna"></a>

1. **Alignment**: AASRA tool was used for alignment.
2. **Counting**: FeatureCounts with default parameters using saf file produced by AASRA.
3. **Reference**: Known Rattus novergicus miRNAs from miRBase (release 22.1), piRNAs from piRbase, tsRNAs from GtRNAdb/rn6, and rsRNAs from Ensembl (ncRNAs, Rnor_6.0 release 104).
4. **Differential Expression Analysis**: DESeq2 used with sncRNAs having at least 20 counts. Adjusted p-value < 0.1 and fold change > 2 were thresholds.

## Scripts and Reference Data

The scripts for lncRNA-Seq analysis:
- `code_lncRNA.sh` #Bash script with the command lines of initial steps of pipeline (Quality Control,Alignment,Quantification).
- `lncRNA_DEA.Rmd` #R script for filtering and Differential Expression Analysis.
- `lncRNA_Reference_info.txt` #Text file with information about reference used for lncRNA analysis.

The scripts for sncRNA-Seq analysis:
- `code_sncRNA.sh` #Bash script with the command lines of initial steps of pipeline (Quality Control,Alignment,Quantification).
- `sncRNA_DEA.Rmd` #R script for filtering and Differential Expression Analysis.
- `sncRNA_Reference_info.txt` #Text file with information about reference used for sncRNA analysis.



## Dependencies

- Hisat2
- FeatureCounts
- DESeq2
- AASRA
- R (with relevant libraries)
- CPAT

## Usage

1. Install proper all dependencies.
2. Follow the provided bash and R scripts for each analysis using proper directories and references.
3. Ensure all dependencies are installed and properly configured and download the specie specific references from ensembl, gencode or specific small RNAs databanks.

