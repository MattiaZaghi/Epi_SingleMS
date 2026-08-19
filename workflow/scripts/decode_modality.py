import argparse
import gzip
from Bio import SeqIO

def is_valid_sequence(seq):
    valid_bases = set('ACGT')
    return set(seq).issubset(valid_bases)

def demultiplex_fastq(input_file, barcode, output_file):
    with gzip.open(input_file, 'rt') as input_handle, gzip.open(output_file, 'wt') as output_handle:
        for record in SeqIO.parse(input_handle, 'fastq'):
            read = str(record.seq)
            qual = record.letter_annotations["phred_quality"]
            barcode_seq = read[0:8]
            if barcode_seq == barcode and is_valid_sequence(read):
                new_record = record[:]
                new_record.letter_annotations["phred_quality"] = qual
                SeqIO.write(new_record, output_handle, 'fastq')

def main():
    # Create command-line argument parser
    parser = argparse.ArgumentParser(description='Demultiplex FASTQ files based on barcode.')
    parser.add_argument('-iR2', '--input1', dest='input_file1', required=True, help='First input FASTQ file (.gz)')
    parser.add_argument('-iR1', '--input2', dest='input_file2', required=True, help='Second input FASTQ file (.gz)')
    parser.add_argument('-b', '--barcode', dest='barcode', required=True, help='Barcode sequence')
    parser.add_argument('-oR2', '--output1', dest='output_file1', required=True, help='First output FASTQ file (.gz)')
    parser.add_argument('-oR1', '--output2', dest='output_file2', required=True, help='Second output FASTQ file (.gz)')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call demultiplex_fastq function for input1 and output1
    demultiplex_fastq(args.input_file1, args.barcode, args.output_file1)

    # Process input2 based on output1 and write to output2
    index_map = {}
    with gzip.open(args.output_file1, "rt") as file:
        for record in SeqIO.parse(file, "fastq"):
            header = record.id.rsplit("/", 1)[0]  # strip /2
            index_map[header] = True

    with gzip.open(args.input_file2, "rt") as input_file, gzip.open(args.output_file2, "wt") as output_file:
        for record in SeqIO.parse(input_file, "fastq"):
            header = record.id.rsplit("/", 1)[0]  # strip /1
            if header in index_map:
                SeqIO.write(record, output_file, "fastq")

if __name__ == "__main__":
    main()

