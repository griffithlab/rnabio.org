---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Strand Settings
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-12-01
---

### Strand Settings

| Property / Tool                                      | RF Stranded                        | FR Stranded                     | Unstranded                 |
|------------------------------------------------------|------------------------------------|---------------------------------|----------------------------|
| **Library Kit Examples**                             | TruSeq Strand Specific Total RNA   | NuGEN Encore                    | NuGEN OvationV2            |
| **Stranded?**                                        | Yes                                | Yes                             | No                         |
| **5p to 3p IGV**                                     | F2R1                               | F1R2                            | F2R1 or F1R2               |
| **TopHat (--library-type parameter)**                | fr-firststrand                     | fr-secondstrand                 | fr-unstranded              |
| **HISAT2 (--rna-strandness)**                        | R/RF                               | F/FR                            | NONE                       |
| **HTSeq (--stranded/-s)**                            | reverse                            | yes                             | no                         |
| **Picard (CollectRnaSeqMetrics STRAND_SPECIFICITY)** | SECOND_READ_TRANSCRIPTION_STRAND   | FIRST_READ_TRANSCRIPTION_STRAND | NONE                       |
| **Kallisto quant**                                   | --rf-stranded                      | --fr-stranded                   | NONE                       |
| **StringTie**                                        | --rf                               | --fr                            | NONE                       |

