This document details the uses of each script in this directory.
Many scripts refer to directories and files with names containing FM and P.
FM = Functionally Monogynous (OT3b assembly); P = Polygynous (PF assembly).
NOTE: swap directories and files out as needed.
For any questions email jpstevenson2018@outlook.com

------------------------------------------------------------------------------

##BLAST_aa_FM.script 
BLASTp script blasting amino acid sequence data from BRAKER output against \
swissprot database. (To use the same script for P data, replace the -query)

##BLAST_codingseq_FM.script
BLASTn script blasting nucleotide sequence data from BRAKER output against \
NCBI nt database. (To use the same script for P data, replace the -query)

##BRAKER_FM.script ; BRAKER_P.script
Script used for BRAKER2 with FM data, internally calling the Anaconda \
environment. Input files (--genome; --bam) can be replaced.
NOTE: Might want to add in --alternatives_from_evidence command.

##clean_results_genes_FM.script ; clean_results_genes_P.script
Script calling getAnnoFastaFromJoingenes.py which finds predicted genes in \
Augutus FM output with internal stop codons producing a list of these genes \
and codingseq and aa output files with these bad genes removed.

##L_acer_data_FASTQC.script 
FASTQC script used on RNA-seq data.

##L_acer_Trim_default_quality.script
TrimGalore script used on RNA-seq data.
  
##SCRMshaw_FM.script
SCRMshaw script for predicting CRMs from Augustus output gff3 file. Script \
internally calls Anaconda environment. Uses training data whose directories \
should be listed in a separate file and fed in with --traindirlst. (To use \
the same script for P data, replace the --gff)

##SCRMshaw_top_hits.script
Script to find the 250 top hits from SCRMshaw results.

##sort_merge.sh
Simple bash script for sorting and merging the top hits.

##STAR_alignment_FM.script ; STAR_alignment_P.script
Script for aligning RNA-seq data with genome assembly index from the \
Genome_Generate script.

##STAR_Genome_Generate_FM.script ; STAR_Genome_Generate_P.script
Script for generating an index of your genome assembly used when aligning.

##Trinity_L_acer_FM.script ; Trinity_L_acer_P.script
Script for genome guided transcriptome assembly using Trinity and the STAR \
alignment file.
