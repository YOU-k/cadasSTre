## **Overview**
cadasSTre is a cross-platform data for sequencing-based ST benchmarking which allows for the evaluation of the performance of each technology in terms of spatial resolution, capture efficiency, and molecular diffusion. 

This is an overview of a set of scripts used for processing various types of data, including their execution order, input & output files, and a brief description for each script.

## **How do I use these scripts for my own analysis?**

The analysis pipeline follows a certain script execution order. To process the data, it is recommended to execute these scripts in the following order:
1. Begin with bam_barcode_matching.py to clean the BAM file, perform bead barcode matching, and generate a whitelist
2. Next, run dbit.R to process genome alignment data, perform barcode detection, and counting
3. Use visium_eb1.sh to process 10x Genomics Visium data
4. Employ stereoseq.py to process other genome data and generate a gene-UMI matrix
5. Finally, execute select_reads_copy.py to extract gene tags, UMI tags, and generate a gene-UMI matrix from a BAM file


## **What do these scripts do and their corresponding input and output files?**

Script 1: bam_barcode_matching.py
Performs BAM file cleanup, bead barcode matching, and generates a whitelist based on barcode information.
Input: a BAM file and a reference data file defining barcodes
Output: cleaned BAM file and the generated whitelist

Script 2: dbit.R
Used for processing genome alignment data, barcode detection, and counting.
Input: sample ID, data file path, reference genome file, annotation file, and an index directory
Output: gene count information, annotation file, and a comparison file

Script 3: visium_eb1.sh
Processes 10x Genomics Visium data using the Space Ranger tool.
Input: sample ID, transcriptome data, FASTQ file path, image file, slide and area information, local thread, and memory settings
Output: results of processing 10x Genomics Visium data

Script 4: stereoseq.py
Extracts gene tags, UMI tags, and other relevant information from a BAM file and generates a gene-UMI matrix.
Input: BAM file, gene tags, UMI tags, and sub-property definitions
Output: feature file, barcode file, and matrix file

Script 5: select_reads_copy.py
Extracts gene tags, UMI tags, and other relevant information from a BAM file and generates a gene-UMI matrix.
Input: BAM file, gene tags, UUMI tags, and sub-property definitions
Output: feature file, barcode file, and matrix file

The paths and names of input and output files can be customized as needed. These scripts are designed to help you process different types of data, from cleaning and matching bead barcodes to generating gene-UMI matrices, providing convenient data preparation tools for subsequent analysis. If you need more information or have any questions, please feel free to reach out.
