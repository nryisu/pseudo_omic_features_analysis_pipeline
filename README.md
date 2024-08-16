# pseudo_omic_features_analysis_pipeline
 This repository processes the pseudo-omic features induced by the extended S phase in cell composition.
 
# Getting started
The workflow is described in the Figure folder.  
   The Figure 1 section corresponds to the methods of WGS-related processes, such as the pipeline of WGS data, simulation S phase ratio, CNV segment calling, and visualization...  
   The Figure 2 section corresponds to the methods of NicE-seq data analysis.  
   The Figure 3 section corresponds to the methods of the DNA methylation data analysis.  
   The Figure 4 section corresponds to the methods of the single-cell RNA-seq analysis.  
   For details, see the documentation of each of the scripts below the folder.  

# Dependencies
You will need the following software:
     TrimGalore v0.6.3, BWA v0.7.17, DNAcopy v1.58.0, STAR v2.7.1a, featureCounts v1.6.4, MAST v1.10.0,  MACS2 v2.1.1, sva v3.35.2

# References
 1. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 2013, 29(1): 15-21.
 2. Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 2013, 30(7): 923-930.
 3. Krueger F. Trim Galore: a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation  Bisufite-Seq) libraries. URL http://www bioinformatics babraham ac uk/projects/trim_galore/(Date of access: 28/04/2016) 2012.
 4. McCarthy DJ, Campbell KR, Lun AT, Wills QF. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 2017, 33(8): 1179-1186.
 5. McDavid A, Finak G, Yajima M. MAST: model-based analysis of single cell transcriptomics. R package version 0931 2015.
 6. Seshan VE, Olshen AB. DNAcopy: a package for analyzing DNA copy data. Bioconductor Vignette 2014: 1-7.
 7. Weddington N, Stuy A, Hiratani I, Ryba T, Yokochi T, Gilbert DM. ReplicationDomain: a visualization tool and comparative database for genome-wide replication timing data. BMC bioinformatics 2008, 9(1): 530.
