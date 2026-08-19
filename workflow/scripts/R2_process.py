
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-iR2", "--input", required=True, help="input file")
ap.add_argument("-oR1", "--output_R1", required=True, help="output file R1")
ap.add_argument("-oR2", "--output_R2", required=True, help="output file R2")
args = vars(ap.parse_args())

input_file_R1 = args["input"]
output_file_R1 = args["output_R1"]
output_file_R2 = args["output_R2"]

#seq_start = 125  # 22bp primer + 8bp BC2 + 30bp linker2 + 8bp BC1 + 30bp linker1 + 8bp modality barcode + 19bp ME (chemV2 with UMI)
seq_start = 8  # 8bp modality barcode + 80bp genomic DNA + 8bp adaptor sequences + 16bp single cell barcodes
seq_end = 88  # 

bc_start = 96
bc_end = 112

with gzopen(input_file_R1, "rt") as in_handle_R1, gzopen(output_file_R1, "wt") as out_handle_R1, gzopen(output_file_R2, "wt") as out_handle_R2:
    for title, seq, qual in FastqGeneralIterator(in_handle_R1):
        new_seq_R1 = seq[seq_start:seq_end]
        new_qual_R1 = qual[seq_start:seq_end]

        barcode = seq[bc_start:bc_end]  # !!! last 16nt
        new_qual_R2 = qual[bc_start:bc_end]
        out_handle_R1.write("@%s\n%s\n+\n%s\n" % (title, new_seq_R1, new_qual_R1))
        out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual_R2))

