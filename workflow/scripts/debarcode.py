import os
import argparse
from glob import glob
import sys
from re import split
import regex
import gzip
from contextlib import ExitStack
from collections import defaultdict
import time

import yaml
from pysam import FastxFile
import Levenshtein

def log(message):
    sys.stderr.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + message + "\n")

class bcdCT:
    def __init__(self,args):
        self.platform = args.platform
        if self.platform == 'mgi' and not args.single_cell:
            log("*** Error: --platform mgi requires --single_cell (there is no cell barcode to extract otherwise) ***")
            sys.exit(1)
        self.detect_input(args.input)
        self.detect_reads()
        self.single_cell=args.single_cell
        self.out_prefix=args.out_prefix
        if self.single_cell:
            self.out_reads = ['R1','R2','R3']
        else:
            self.out_reads = ['R1','R3']

        if args.name:
            self.name = args.name
        else:
            self.autodetect_name()


        self.autodetect_barcodes(args)
        self.prep_out_filenames()

    def detect_input(self,input):
        Error_message="*** Error: Wrong input files specified. The input must be either folder with _R1_*.fastq.gz _R2_*.fastq.gz _R3_*.fastq.gz files or paths to the files themselves ***" +\
        "The files should be placed in the same folder" +\
        "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"

        input = [os.path.abspath(x) for x in input]
        if len(input) == 1 and os.path.isdir(input[0]):     # Case input is single directory
            self.input_dir = input[0]
            self.input_files = []
            self.input_files.extend(glob(self.input_dir + "/*.fastq.gz"))
            self.input_files.extend(glob(self.input_dir + "/*.fq.gz"))

        elif len(input) > 1:                                  # Case input are multiple files
            self.input_files = input
            self.input_dir = list(set([os.path.dirname(x) for x in self.input_files]))
            if not len(self.input_dir) == 1:
                log(Error_message)
                sys.exit(1)
            if not sum([x.endswith('.fastq.gz') or x.endswith('.fq.gz') for x in self.input_files]) == len(self.input_files):
                sys.exit(1)
                log(Error_message)
        else:
            sys.exit(1)
            log(Error_message)

    def detect_reads(self):
        Error_message="*** Error: Please specify exactly one _R1_ _R2_ and _R3_ file or folder with exactly one of each files ***" + \
                      "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"
        Error_message_mgi="*** Error: Please specify exactly one _R1_ and _R2_ file (MGI platform, no _R3_) or folder with exactly one of each files ***" + \
                      "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/"
        self.path_in = {}
        self.path_in['R1'] = [x for x in self.input_files if "_R1_" in x]
        self.path_in['R2'] = [x for x in self.input_files if "_R2_" in x]

        if self.platform == 'mgi':
            if len(self.path_in['R1']) != 1 or len(self.path_in['R2']) != 1:
                log(Error_message_mgi)
                sys.exit(1)
        else:
            self.path_in['R3'] = [x for x in self.input_files if "_R3_" in x]
            if len(self.path_in['R1']) != 1 or len(self.path_in['R2']) != 1 or len(self.path_in['R3']) != 1:
                log(Error_message)
                sys.exit(1)

        self.path_in = {key:self.path_in[key][0] for key in self.path_in.keys()}

        if self.platform == 'mgi':
            # No physical R3 file on MGI runs; R3 output naming is derived from R2's
            # (prep_out_filenames() needs a path_in['R3'] entry for every out_read).
            self.path_in['R3'] = self.path_in['R2'].replace('_R2_', '_R3_')

    def in_handles(self,stack):
        reads = ['R1','R2'] if self.platform == 'mgi' else ['R1','R2','R3']
        in_stack = {x: stack.enter_context(FastxFile(self.path_in[x],'r')) for x in reads}
        return in_stack

    def prep_out_filenames(self):
        self.path_out = {barcode: {} for barcode in self.picked_barcodes}
        self.path_out  = {barcode: {read: "{0}/barcode_{1}/{2}".format(self.out_prefix,barcode,os.path.basename(self.path_in[read])) for read in self.out_reads} for barcode in self.picked_barcodes}
        # If args.name is specified, replace the sample_id prefix with the one specified in args.name
        # e.g. nanoCT_MB22_001_S1_L001_R1_001.fastq.gz is input --name is test
        # Change to test_S1_L001_R1_001.fastq.gz
        if args.name:
            for barcode in self.path_out:
                for read in self.path_out[barcode]:
                    sample_id = split('_S[0-9]+_', os.path.basename(self.path_out[barcode][read]))[0].strip("_")
                    self.path_out[barcode][read] = self.path_out[barcode][read].replace(sample_id,args.name)

    def autodetect_name(self):
        Error_message = "*** Error: Prefix for R1 R2 R3 files not the same. Please use the same prefix for all the files or specify experiment name ***"

        self.name = [split("_R[0-9]_", str(x)) for x in self.path_in.values()]
        self.name = [x[0] for x in self.name]

        if len(list(set(self.name))) > 1:
            log(Error_message)
            sys.exit(1)

        self.name = self.name[0].split("/")[-1]

    def create_out_handles(self,stack):
        for bcd in self.picked_barcodes:
            os.makedirs(self.out_prefix + "/barcode_" + bcd, exist_ok=True)
        self.out_stack = {barcode: {read: stack.enter_context(gzip.open(self.path_out[barcode][read],'wt'))for read in self.out_reads} for barcode in self.picked_barcodes}


    def __iter__(self):
        if self.platform == 'mgi':
            with FastxFile(self.path_in['R1']) as f1, FastxFile(self.path_in['R2']) as f2:
                for r1,r2 in zip(f1,f2):
                    yield r1, r2
        else:
            with FastxFile(self.path_in['R1']) as f1, FastxFile(self.path_in['R2']) as f2, FastxFile(self.path_in['R3']) as f3:
                for r1,r2,r3 in zip(f1,f2,f3):
                    yield r1, r2, r3

    def autodetect_barcodes(self,args):
        barcodes = defaultdict(int)
        n=0
        if self.platform == 'mgi':
            for read1,read2 in self:
                # MGI: modality barcode, insert and cell barcode are fused on R2
                # at fixed offsets (modality barcode = first barcode_length nt).
                read_barcode = parse_mgi_read(
                    read2.sequence, args.barcode_length, args.mgi_insert_length,
                    args.mgi_adaptor_length, args.cell_barcode_length)[0]
                if read_barcode is None:
                    continue
                barcodes[read_barcode] += 1
                n += 1
                if n == 50000:
                    break
        else:
            for read1,read2,read3 in self:
                # V4: the modality barcode is on R3 (sequencing R2), flanked by two MEs.
                read_barcode, trim_end, orientation = find_modality_barcode(
                    read3, args.me, args.barcode_length, args.me_mismatch)
                if read_barcode is None:
                    continue
                barcodes[read_barcode] += 1
                n += 1
                if n == 50000:
                    break

        top_barcodes = sorted(barcodes, key=barcodes.get, reverse=True)[:args.Nbarcodes]
        picked_barcodes = {key: barcodes[key] for key in top_barcodes}
        log("Detected following most abundant barcodes out of first {} barcodes:\n{}".format(n, picked_barcodes))
        if args.barcode != "None":
            self.picked_barcodes = args.barcode
            log("Barcode specified for demultiplexing [{barcode}] in top found barcodes: {bool} ".format(bool = [(x,x in picked_barcodes.keys()) for x in args.barcode], barcode = args.barcode))
        else:
            self.picked_barcodes = [i for i in picked_barcodes.keys()]

        if args.report_no_hit:
            self.picked_barcodes.append('no_ME')
            log("Reads with no ME-flanked modality barcode will be reported (bucket 'no_ME') due to --report_no_hit flag")
        
        print('final barcodes used for demultiplexing:')
        print(self.picked_barcodes)
        

def find_modality_barcode(read, me, barcode_length, me_mismatch):
    """Locate the ME-flanked modality barcode on a read (V4 library design).

    R3 (sequencing R2) layout, 5'->3':
        [handle] - ME - <barcode_length nt modality barcode> - ME - <genomic insert>

    Returns (barcode, trim_end, orientation) where ``trim_end`` is the index just
    past the SECOND ME (the start of the genomic insert) in the strand of the
    returned orientation, and ``orientation`` is '+' (as sequenced) or '-'
    (reverse complement). Returns (None, None, None) when no valid cassette is
    found. The barcode is located by anchoring on the ME motif on BOTH sides
    rather than by a fixed offset, so a variable read-start position is
    tolerated. If the forward strand yields no cassette the reverse complement
    is tried, so orientation is read from the data instead of assumed.
    """
    barcode, trim_end = locate_cassette(read.sequence, me, barcode_length, me_mismatch)
    if barcode is not None:
        return barcode, trim_end, '+'
    barcode, trim_end = locate_cassette(revcompl(read.sequence), me, barcode_length, me_mismatch)
    if barcode is not None:
        return barcode, trim_end, '-'
    return None, None, None


def locate_cassette(seq, me, barcode_length, me_mismatch):
    """Find  ME - <barcode> - ME  on a single strand. Returns (barcode, trim_end)."""
    first = find_seq_first(me, seq, nmismatch=me_mismatch)              # left flanking ME
    if first is None:
        return None, None
    bc_start = first + len(me)
    bc_end   = bc_start + barcode_length
    barcode  = seq[bc_start:bc_end]
    if len(barcode) < barcode_length:                                  # read too short for a full barcode
        return None, None
    pos = find_seq_first(me, seq[bc_end:], nmismatch=me_mismatch)       # right flanking ME
    if pos is None or pos > me_mismatch:                               # must sit immediately after the barcode
        return None, None
    trim_end = bc_end + pos + len(me)                                  # start of the genomic insert
    return barcode, trim_end


def parse_mgi_read(seq, barcode_length, insert_length, adaptor_length, cell_barcode_length):
    """Split a fused MGI R2 read into modality barcode, genomic insert and cell
    barcode by FIXED offsets (no mosaic-end anchoring -- MGI reads are fixed
    length and carry no ME).

    R2 layout, 5'->3' (default lengths shown):
        [0 : bc]                       modality barcode      (bc = 8)
        [bc : bc+insert]               genomic insert        (insert = 80)  -> [8:88]
        [bc+insert : +adaptor]         adaptor, discarded    (adaptor = 8)  -> [88:96]
        [cell_start : cell_start+cbc]  cell barcode          (cbc = 16)     -> [96:112]

    Mirrors the reference R2_process.py / decode_modality.py offsets. The modality
    barcode is the bare first ``barcode_length`` nt; matching to the known set
    (with mismatch tolerance) is done downstream by ``match_barcode``.

    Returns (modality_barcode, cell_start, cell_end, insert_start, insert_end),
    or a 5-tuple of None if the read is too short to hold a full cell barcode.
    """
    insert_start = barcode_length
    insert_end   = barcode_length + insert_length
    cell_start   = insert_end + adaptor_length
    cell_end     = cell_start + cell_barcode_length
    if len(seq) < cell_end:                                            # read too short for a full cell barcode
        return None, None, None, None, None
    modality_barcode = seq[0:barcode_length]
    return modality_barcode, cell_start, cell_end, insert_start, insert_end


def match_barcode(read_barcode, picked_barcodes, mismatch, statistics):
    """Match ``read_barcode`` against ``picked_barcodes`` within ``mismatch``
    edits, bumping the same statistics counters the illumina path always has.
    Returns the matched barcode, or None on no match / ambiguous match."""
    read_barcode_distance = {barcode: Levenshtein.distance(read_barcode,barcode) for barcode in picked_barcodes}
    n_hits = sum([x <= int(mismatch) for x in read_barcode_distance.values()])
    if n_hits == 0:
        # ME cassette found but no barcode match
        statistics["no_barcode_match"] += 1
        return None
    if n_hits > 1:
        # ME cassette found but multiple barcode matches
        statistics["multiple_barcode_matches"] += 1
        return None
    return min(read_barcode_distance,key=read_barcode_distance.get)


def reconstitute_read(read, trim_end, orientation):
    """Rebuild R3 as it would be WITHOUT the modality-barcode cassette: strip
    everything up to and including the second ME, leaving only the genomic
    insert. Mirrors the old R2 cell-barcode trimming convention (reassign
    .sequence/.quality on the record in place, then write str(read))."""
    if orientation == '-':
        read.sequence = revcompl(read.sequence)
        read.quality  = read.quality[::-1] if read.quality else read.quality
    read.sequence = read.sequence[trim_end:]
    read.quality  = read.quality[trim_end:] if read.quality else read.quality
    return read

def revcompl(seq):
    revcomp_table = {
        "A": "T",
        "G": "C",
        "C": "G",
        "T": "A",
        "N": "N"
    }
    complement = "".join([revcomp_table[letter] for letter in seq.upper()])  # Complement
    return complement[::-1]  # Reverse

def rev(seq):
    return seq[::-1]

def strip_mate_suffix(name):
    """Drop a trailing '/1' or '/2' mate suffix from an MGI read name so the
    R1/R2/R3 outputs share an identical name (cellranger pairs reads by name)."""
    if name and len(name) > 2 and name[-2:] in ('/1', '/2'):
        return name[:-2]
    return name

def find_seq(pattern, DNA_string, nmismatch=2):
    for n in range(0,nmismatch + 1):
        r = regex.compile('({0}){{e<={1}}}'.format(pattern, n))
        res = r.finditer(DNA_string)
        hit = [x.start() for x in res]
        if len(hit) == 0:
            continue
        if len(hit) > 1:
            return None
        if len(hit) == 1:
            return int(hit[0])
    return None

def find_seq_first(pattern, DNA_string, nmismatch=2):
    """Return the start of the earliest fuzzy match of ``pattern`` (fewest edits
    first), or None. Unlike ``find_seq`` it does not require the match to be
    unique -- used to anchor on the ME motif, which can recur inside the insert."""
    for n in range(0, nmismatch + 1):
        m = regex.compile('({0}){{e<={1}}}'.format(pattern, n)).search(DNA_string)
        if m:
            return m.start()
    return None

def main(args):
    exp = bcdCT(args)
    statistics = defaultdict(int)
    log("Creating file output handles ")
    with ExitStack() as stack:
        exp.create_out_handles(stack)
        n = 0
        log("Starting demultiplexing ")

        if exp.platform == 'mgi':
            for read1,read2 in exp:
                n+=1
                if n % 5000000 == 0:
                    log("{} reads processed".format(n))
                # MGI names carry /1 /2 mate suffixes; strip so R1/R2/R3 names match for cellranger.
                read1.name = strip_mate_suffix(read1.name)
                read2.name = strip_mate_suffix(read2.name)
                assert (read1.name == read2.name)                                                            # Make sure the fastq files are ok

                # MGI design: modality barcode, genomic insert and cell barcode
                # are all fused on R2 at fixed offsets -- there is no separate R3.
                modality_barcode, cell_start, cell_end, insert_start, insert_end = parse_mgi_read(
                    read2.sequence, args.barcode_length, args.mgi_insert_length,
                    args.mgi_adaptor_length, args.cell_barcode_length)

                if modality_barcode is None:
                    # R2 too short to hold a full cell barcode
                    statistics["too_short_read"] += 1
                    continue

                hit_barcode = match_barcode(modality_barcode, exp.picked_barcodes, args.mismatch, statistics)
                if hit_barcode is None:
                    if args.report_no_hit:
                        hit_barcode = 'no_ME'
                    else:
                        continue

                cell_seq    = read2.sequence[cell_start:cell_end]
                cell_qual   = read2.quality[cell_start:cell_end] if read2.quality else read2.quality
                insert_seq  = read2.sequence[insert_start:insert_end]
                insert_qual = read2.quality[insert_start:insert_end] if read2.quality else read2.quality

                statistics[hit_barcode] += 1
                if hit_barcode in exp.picked_barcodes:
                    # Write the outputs
                    exp.out_stack[hit_barcode]['R1'].write('{}\n'.format(str(read1)))

                    read2.sequence, read2.quality = insert_seq, insert_qual
                    exp.out_stack[hit_barcode]['R3'].write('{}\n'.format(str(read2)))

                    read2.sequence, read2.quality = cell_seq, cell_qual
                    exp.out_stack[hit_barcode]['R2'].write('{}\n'.format(str(read2)))

        else:
            for read1,read2,read3 in exp:
                n+=1
                if n % 5000000 == 0:
                    log("{} reads processed".format(n))
                assert (read1.name == read2.name == read3.name)                                                 # Make sure the fastq files are ok

                # V4 design: the 8-nt modality barcode lives on R3 (sequencing R2),
                # flanked by two Tn5 mosaic ends (ME). The 16-nt cell barcode stays
                # untouched on R2 (i5 index).
                read_barcode, trim_end, orientation = find_modality_barcode(
                    read3, args.me, args.barcode_length, args.me_mismatch)

                if read_barcode is None:
                    # No ME-flanked cassette found on R3
                    statistics["no_ME_cassette"] += 1
                    if args.report_no_hit:
                        hit_barcode = 'no_ME'
                    else:
                        continue
                else:
                    hit_barcode = match_barcode(read_barcode, exp.picked_barcodes, args.mismatch, statistics)
                    if hit_barcode is None:
                        continue

                    # Reconstitute R3 as it would be WITHOUT the modality-barcode cassette
                    # (strip through the second ME, leaving only the genomic insert).
                    read3 = reconstitute_read(read3, trim_end, orientation)

                # The cell barcode is carried on R2 (i5) unchanged; guard truncated reads.
                if exp.single_cell and len(read2.sequence) < args.cell_barcode_length:
                        statistics["too_short_read"] += 1
                        continue

                statistics[hit_barcode] += 1
                if hit_barcode in exp.picked_barcodes:
                    # Write the outputs
                    exp.out_stack[hit_barcode]['R1'].write('{}\n'.format(str(read1)))
                    exp.out_stack[hit_barcode]['R3'].write('{}\n'.format(str(read3)))
                    if args.single_cell:
                        exp.out_stack[hit_barcode]['R2'].write('{}\n'.format(str(read2)))


    # Write the statistics file
    with open("{0}/{1}_statistics.yaml".format(exp.out_prefix,exp.name), 'w') as f:
        yaml.dump(statistics, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="DESCRIPTION: \n\nThis script demultiplexes nano-CUT&Tag (V4 library design) sequencing data by extracting and matching the antibody/modality barcode from R3 (which is the sequencing R2 / genomic read).\nThe script supports both bulk and single-cell data and writes sorted reads into separate output files for each detected barcode.\n""" +
                                     """V4 read structure -- the 8-nt modality barcode is flanked by two Tn5 mosaic ends (ME) on R3:\n""" +
                                     """  R3:  [handle ...TGTCTCGGGGCTCGGAG] - ME - <8nt modality barcode> - ME - <genomic insert>\n""" +
                                     """       ME = AGATGTGTATAAGAGACAG\n""" +
                                     """The barcode is located by anchoring on the ME motif on both sides (not by a fixed offset), and the cassette is trimmed from R3 through the second ME so the aligner receives only the genomic insert.\n""" +
                                     """The 16-nt cell barcode remains on R2 (i5 index) and is passed through unchanged; R1 is passed through unchanged.\n""" +
                                     """(Legacy V2 --pattern/--no_barcode_seq spacer options are still accepted but no longer used.)\n""" +
                                     """--platform mgi mode: for MGI-sequenced libraries only R1+R2 are produced (no R3). R2 is one\n""" +
                                     """fused read carrying the modality barcode, the genomic insert AND the cell barcode at FIXED offsets:\n""" +
                                     """  R2:  <8nt modality barcode> | <80nt genomic insert> | <8nt adaptor> | <16nt cell barcode>   (positions [0:8][8:88][88:96][96:112])\n""" +
                                     """The modality barcode is the bare first 8 nt (no ME); the genomic insert R2[8:88] is written as R3 and the\n""" +
                                     """cell barcode R2[96:112] as R2, so cellranger receives R1(passthrough)/R2(cell barcode)/R3(insert). Offsets are\n""" +
                                     """tunable via --barcode_length/--mgi_insert_length/--mgi_adaptor_length/--cell_barcode_length.\n""" +
                                     """--platform mgi requires --single_cell and only accepts one _R1_ and one _R2_ input file (no _R3_).\n""" +
                                     """ Note: If demultiplexing multiple lanes, run for each lane separetely and then merge the output files before or after alignment""",
                                     usage=""
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC                      # One specific barcode from single-cell data "
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --barcode ATAGAGGC TATAGCCT             # Two specific barcodes from single-cell data "
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --single_cell --Nbarcodes 3                           # Top 3 barcodes from single-cell data without specifying the barcodes - use carefuly and double check"
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz /path/to/input_R3.fastq.gz -o /path/to/output_folder --Nbarcodes 3                                         # Top 3 barcodes from bulk data "
                                           "python debarcode.py -i /path/to/input_R1.fastq.gz /path/to/input_R2.fastq.gz -o /path/to/output_folder --platform mgi --single_cell --barcode ATAGAGGC                                    # MGI platform (R1+R2 only) single-cell data ",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        nargs='+',
                        help='path to input R1,R2,R3 .fastq.gz files [3 files required]')

    parser.add_argument('-o', '--out_prefix',
                        type=str,
                        required=True,
                        help='Prefix to where to put the output files; Diretory will be created')

    parser.add_argument('--me',
                        type=str,
                        default="AGATGTGTATAAGAGACAG",
                        help='Tn5 mosaic-end (ME) motif flanking the modality barcode on R3 (Default: %(default)s)')

    parser.add_argument('--barcode_length',
                        type=int,
                        default=8,
                        help='Length of the modality/antibody barcode in nt (Default: %(default)s)')

    parser.add_argument('--me_mismatch',
                        type=int,
                        default=2,
                        help='Maximum edits allowed when locating each ME motif on R3 (Default: %(default)s)')

    parser.add_argument('-p', '--pattern',
                        type=str,
                        default="CAGACGCG",
                        help='[LEGACY / V2, no longer used] spacer that followed the antibody barcode on R2 \n \
                                  (Default: %(default)s)')

    parser.add_argument('--single_cell',
                        default=False,
                        action='store_true',
                        help='Data is single cell CUT&Tag (Default: %(default)s)')

    parser.add_argument('--name',
                        type=str,
                        default=None,
                        help='Custom name for the experiment (Default: Autodetect from filename)')

    parser.add_argument('--mismatch',
                        type=int,
                        default=1,
                        help='Maximum edit distance allowed when matching the modality barcode to the known set (Default: %(default)s)')

    parser.add_argument('--Nbarcodes',
                        type=int,
                        default=10,
                        help='Number of barcodes in experiment (Default: %(default)s)')

    parser.add_argument('--barcode',
                        type=str,
                        nargs="+",
                        default='None',
                        help='Specific barcode to be extracted [e.g. ATAGAGGC] (Default: All barcodes [see --Nbarcodes])')

    parser.add_argument('--no_barcode_seq', type=str,
                        default='GTGTAGATCTCGGTGGTCGCCGTATCATT',
                        help='[LEGACY / V2, no longer used] sequence indicating unbarcoded reads')

    parser.add_argument('--report_MeA',
                        action='store_true',
                        help='[LEGACY / V2, no longer used] include reads with no barcode and standard MeA sequence in the output')

    parser.add_argument('--report_no_hit',
                        action='store_true',
                        help="Include reads with no ME-flanked modality barcode in the output (bucket 'no_ME')")

    parser.add_argument('--platform',
                        type=str,
                        choices=['illumina','mgi'],
                        default='illumina',
                        help='Sequencing platform / read layout (Default: %(default)s). \n' + \
                             '"mgi" reads only R1+R2 (no R3): R2 carries the modality-barcode\n' + \
                             'cassette, the genomic insert and the cell barcode all fused together.\n' + \
                             'Requires --single_cell.')

    parser.add_argument('--cell_barcode_length',
                        type=int,
                        default=16,
                        help='Length of the single-cell cell barcode in nt (Default: %(default)s)')

    parser.add_argument('--mgi_insert_length',
                        type=int,
                        default=80,
                        help='[--platform mgi] length in nt of the genomic insert on the fused R2 read, \n \
                                  immediately after the modality barcode (Default: %(default)s)')

    parser.add_argument('--mgi_adaptor_length',
                        type=int,
                        default=8,
                        help='[--platform mgi] length in nt of the adaptor between the genomic insert \n \
                                  and the cell barcode on the fused R2 read (Default: %(default)s)')


    args = parser.parse_args()
    log("Starting debarcode.py script ")
    log("Input files: \n{}".format("".join(["    " + i + "\n" for i in args.input])))
    if args.barcode != "None":
        log("Provided barcodes to demultiplex: \n{}".format(args.barcode))
    log("Output prefix: {}/".format(args.out_prefix))
    
    main(args)
