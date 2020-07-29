# MHap_DMR
Include the codes to construct methylation haplotypes and identify DMRs based on methylation haplotypes, respectively.

Requirement:

1. Prepare the genome (such as the hf38.fa), cpgIsland information file

2. Samples are sorted mapped BS-seq data (sam file)



Step1: Running Generate_candidate_regions.pl to generate candidate regions.

     Command: perl Generate_candidate_regions.pl genome.fa cpgIsland.txt
     
     Note: genome.fa can be download from USCS, and cpgIsland information also can be extracted from UCSC.
     
     The format of cpgIsland.txt can be this: id chromosome start end. For example: 585	chr1	10468	11240	
               
Step2: Running Haplotype_CpG_Cluster.pl to construct methylation haplotypes for a sample on candidate regions
     
     Command: perl Haplotype_CpG_Cluster.pl sorted_sam_file, CpG_cluster_file
     
     Note: sorted_sam_file is generated by Bismark and sorted by samtools; CpG_cluster_file is a output file of Step 1 after running Generate_candidate_regions.pl.
      
Step3: Running Haplotype_DMR_matrix.pl to identify DMRs.
     
     Command: perl Haplotype_DMR_matrix.pl list.file
     
     Note: list.file is a file containing the filenames of the methylation haplotype files of samples generated by Step2. This step will identify DMRs for the sampples in the list.file.
