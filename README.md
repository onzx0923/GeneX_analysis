# GeneX Variant and Regulatory Network Analysis

This repository contains a comprehensive three-phase bioinformatics workflow developed to investigate the molecular function and regulatory consequences of mutation in Gene X, a putative epigenomic factor. The project integrates genomic, epigenomic and transcriptomic analyses to characterize both the direct physical targets and downstream transcriptional effects associated with disruption of Gene X activity.  

Phase 1 focuses on sequencing quality control and genomic variant identification. Raw paired-end Illumina FASTQ reads were assessed using FastQC, aligned to the Gene X reference sequence using BWA-MEM and processed with SAMtools. Variants were subsequently identified using FreeBayes under a haploid model, followed by filtering and consensus sequence construction using BCFtools.  

Phase 2 involves epigenomic target identification through CUT&RUN analysis. Gene annotation features from the reference GFF3 file were intersected with Gene X CUT&RUN enrichment coordinates using BEDTools to identify genes physically associated with Gene X binding sites.  

Phase 3 integrates the genomic and regulatory datasets to characterize the downstream functional consequences of Gene X mutation. By combining mutation analysis, physical DNA-binding informatio, and transcriptomic data, this workflow enables investigation of the regulatory role of Gene X and its contribution to altered gene expression networks.
