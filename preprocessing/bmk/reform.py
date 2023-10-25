### This script extracts paired sequences containing specific linker sequences from two FASTA files and write them into teo new FASTQ files ###
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

### input file path ###
file2 = '34-1_1.clean.fq.gz'
file1 = '34-1_2.clean.fq.gz'

### output file paths for matching and non-matching pairs ###
output_match2 = 'match_fastq1.fastq.gz'
output_match1 = 'match_fastq2.fastq.gz'
output_no_match2 = 'no_match_fastq1.fastq.gz'
output_no_match1 = 'no_match_fastq2.fastq.gz'

### define linker sequences and possible starting positions ###
linker1 = "TACGT"
linker2 = "CGACTC"
linker1_start_positions = [17,18,19,20,21]  # possible starting positions of linker1
umi_len = 12

### initialize counter ###
counter = 0

### gzip.open handle gzip-compressed files; rt for text mode ### 
with gzip.open(file1, "rt") as handle1, gzip.open(file2, "rt") as handle2, \
         gzip.open(output_match1, "wt") as out_handle_match1, gzip.open(output_match2, "wt") as out_handle_match2, \
         gzip.open(output_no_match1, "wt") as out_handle_no_match1, gzip.open(output_no_match2, "wt") as out_handle_no_match2:
    ### parse fastq from input files ###
    records1 = SeqIO.parse(handle1, "fastq")
    records2 = SeqIO.parse(handle2, "fastq")

    ### iterate over each pair of records ###
    for record1, record2 in zip(records1, records2):
        seq = str(record2.seq)
        match_found = False

        ### iterate over all possible starting positions for linker 1 ###
        for linker1_start in linker1_start_positions:
            linker2_start = linker1_start + len(linker1) + 19

            ### if linker1 and 2 are found ###
            if seq[linker1_start: linker1_start+len(linker1)] == linker1 and seq[linker2_start: linker2_start+len(linker2)] == linker2:
                ### extract barcodes and UMI ###
                bc1 = seq[linker1_start-17:linker1_start]
                bc2 = seq[linker1_start+len(linker1): linker2_start]
                bc3 = seq[linker2_start+len(linker2): linker2_start+len(linker2)+19]
                umi = seq[linker2_start+len(linker2)+19: linker2_start+len(linker2)+19+umi_len]
                new_seq = bc1 + bc2 + bc3 + umi
                ### create a new SeqRecord ### 
                new_record = SeqRecord(Seq(new_seq), id=record2.id, description="", letter_annotations={"phred_quality": record2.letter_annotations["phred_quality"][:len(new_seq)]})
                ### writing to matching files ###
                SeqIO.write(new_record, out_handle_match2, "fastq")
                SeqIO.write(record1, out_handle_match1, "fastq")
                match_found = True
                break  # If a match is found, no need to check other starting positions

        ### no match found -> write records to non-matching files ###
        if not match_found:
            SeqIO.write(record2, out_handle_no_match2, "fastq")
            SeqIO.write(record1, out_handle_no_match1, "fastq")
        ### print process every 10000reads ###
        counter += 1
        if counter % 10000 == 0:
            print(f"Processed {counter} reads")
