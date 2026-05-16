# nanopore_workflow1
# Nanopore Direct RNA-seq Analysis Pipeline

## Overview

This repository contains a Nanopore direct RNA sequencing analysis workflow for:

1. **PolyA tail length extraction**
2. **RNA modification site detection using Dorado + modkit**
3. **m6A-specific analysis using nanopolish + m6Anet**

The workflow supports two complementary strategies:

- **Workflow A:** `Dorado -> minimap2 -> BAM -> polyA extraction + modkit`
- **Workflow B:** `nanopolish -> eventalign -> m6Anet`

---

## Analysis Summary

### Workflow A: Dorado / modkit / PolyA
This workflow is used to:

- retain Dorado-generated tags from FASTQ (`pt`, `MM`, `ML`)
- align reads to the reference transcriptome
- extract **polyA tail length**
- detect **multiple RNA modification types** using `modkit`

### Workflow B: nanopolish / m6Anet
This workflow is used for:

- generating signal-level alignment (`eventalign`)
- preparing m6A candidate features
- running **m6A prediction** using the pretrained `m6Anet` model

---

## Directory Structure

```text
project/
├── 01_alignment/
│   ├── align_reads.sh
│   ├── *.bam
│   └── *.bai
├── 02_polya/
│   ├── extract_polya.sh
│   ├── extract_polya_transcript.sh
│   ├── *_polya_lengths.txt
│   └── *_polya_per_transcript.txt
├── 03_modkit/
│   ├── *.bed
│   ├── *.log
│   └── modkit_pileup.log
├── 04_m6anet/
│   ├── *_summary.txt
│   ├── *_eventalign.txt
│   ├── dataprep/
│   └── run/
├── 05_visualization/
│   ├── polya_boxplot_4groups.R
│   ├── polya_boxplot_4groups.py
│   ├── polya_comparison_4groups.pdf
│   └── polya_comparison_4groups.png
└── README.md
