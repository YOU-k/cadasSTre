### This script aims to process alignment, barcode detection, and counting ###

library(scPipe)
library(SingleCellExperiment)
library(Rsubread)

# use command line argyment to input sampleid
args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]

# set up file path
data_dir = file.path("/path/to/your/data/",sampleid)
fa_fn  = "/path/to/your/sample_assembly.fa"
gff3_fn="/path/to/your/GRCm39.109.gff3"
index_dir="/path/to/your/sample/index"

# define read structure 
read_structure <- list(
  bs1 = -1,   # barcode start position in fq_R1, -1 indicates no barcode
  bl1 = 0,    # barcode length in fq_R1, 0 since no barcode present
  bs2 = 0,    # barcode start position in fq_R2
  bl2 = 16,   # barcode length in fq_R2
  us = 16,    # UMI start position in fq_R2
  ul = 10     # UMI length
)

# alignment
print("###   alignment")
align(
  index = file.path(index_dir),
  readfile1 = file.path(data_dir,"combined_matched.fastq"),
  output_file = file.path(data_dir,"out_aln.bam"),
  nthreads = 32
)

# barcode detection
print("###   detect barcode")
sc_detect_bc(
  infq = file.path(data_dir,"combined_matched.fastq"),
  outcsv = file.path(data_dir, "barcode_anno.csv"),       # bacode annotation output file name
  bc_len = read_structure$bl2, 
  max_reads = 5000000, 
  number_of_cells = 5000,
  white_list_file="/path/to/your/whitelist/first_column.txt"
)

# exon mapping and count
print("###   map to exon and count")
sc_count_aligned_bam(
  inbam = file.path(data_dir,"out_aln.bam"),
  outbam = file.path(data_dir,"out_map.bam"),
  annofn = gff3_fn,
  bc_len = read_structure$bl2,
  UMI_len = read_structure$ul,
  outdir = data_dir,
  bc_anno = file.path(data_dir, "barcode_anno.csv")
)