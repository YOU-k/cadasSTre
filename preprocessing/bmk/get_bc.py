from Bio import SeqIO

# extract sequence from each barcode in fasta file
barcode_to_sequence = {}
with open("/path/to/your/data/sample.fa", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        barcode = record.id  # extract barcode
        sequence = str(record.seq)  # extract corresponding sequence
        barcode_to_sequence[barcode] = sequence

# read barcode combination and write into new file
with open("/path/to/your/data/sample.barcode.tsv", "r") as handle, open("/path/to/your/data/sample/bc_out.tsv", "w") as output_handle:
    for line in handle:
        bc_combination = line.strip()  # extract barcode combination
        bc_list = bc_combination.split('-')  # segment into single barcode
        sequences = [barcode_to_sequence[bc] for bc in bc_list]  # search sequence for each barcode
        merged_sequence = ''.join(sequences)  # combine all sequence
        output_line = bc_combination + '\t' + '\t'.join(sequences) + '\t' + merged_sequence  # generate new line, including barcode combination, seuqnece
        output_handle.write(output_line + '\n')  # write into new file
