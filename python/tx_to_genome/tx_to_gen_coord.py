'''
Goal: Translate a (0-based) transcript coordinate to a (0 based) genome coordinate 
given the following information: chromosome, genomic start position, CIGAR string 
(which indicates how the transcript maps to the genome sequence), and the transcript
coordinate to query.

1. Assumptions:
    * Transcript is always mapped from genomic 5' -> 3'
    * Do not need to accomodate rarer CIGAR annotations (S, H, =, etc.)
    * Empty CIGAR strings are not valid input 
    * Negative transcript coords are not allowed
    * Tab-separated input is provided  

2. Strengths:
    * Checks for valid input (number of columns in input files, 
        CIGAR string formation, empty strings)
    * Run time dependent on CIGAR length but function stops when tr_coord of interest is reached
    (so we don't have to go through the whole CIGAR string) 

3. Weaknesses:

    * My solution does not accomodate the edge case where the CIGAR has >2 letters in a row 
    * Line-by-line parsing does not scale well; could use pandas.read_csv, 
    which is much faster than csvreader under the hood. 
    * RegEx adds some additional overhead; could just iterate through the CIGAR string as needed
    to ensure valid formatting
    * Right now, my script does not allow output of valid queries + transcripts, and instead once
    it runs into a malformed combo, the script stops running.

4. Testing:
    * I tested combinations of valid queries with invalid transcript_files and vice versa (see testing/input/*):
        *I check that CHR column starts with 'CHR'
        *I check that genome_position is an integer
        *I check that no fields are empty
        *I check validity of CIGAR string
        *I check that tr_coord is not negative
        *I check for duplicate transcript_ids

5. Improvements
    * Speed: use pandas to read-in data, refactor to use .applymap() with get_genome_pos()
    * Explicit Unit testing: instead of test cases in input files, write explicit unit tests 

'''
from argparse import ArgumentParser
import csv
import re

def assemble_transcript_dict(transcript_file: str) -> dict:
    '''
    Assembles dict from the transcript_file -- keyed by transcript IDs.

    Args:
        transcript_file: path to the .tsv of file transcripts.
            Expecting 4 columns: transcript ID, genomic chromosome, genomic position, and
            CIGAR string.
    Returns:
        transcripts: a dict of dicts of transcripts, chromosome location,
        start position and CIGAR string.
    Raises:
        ValueError, if there are duplicate transcript IDs.
    '''

    transcripts = {}
    with open(transcript_file, 'r', encoding='utf-8') as tr_f:
        for line in tr_f:
            fields = [i.strip() for i in line.split('\t')]
            if len(fields) != 4:
                raise ValueError(f'Unexpected file format, expect 4 columns, \
                                received {len(fields)} columns.')
            tx_id, chrom, pos, cigar = fields
            if chrom == '' or tx_id == '' or pos == '' or cigar == '':
                raise ValueError(f'Some data is missing, provided TX_ID {tx_id}, \
                                    CHROM {chrom}, POS {pos}, CIGAR {cigar}')
            if not chrom.startswith('CHR'):
                raise ValueError(f'Malformed CHROM name, {chrom}, \
                    expect CHR preceding a number/X/Y')
            if not pos.isnumeric():
                raise ValueError('Genomic position is not numeric!')
            if not is_cigar_valid(cigar):
                raise ValueError(f'Invalid CIGAR string {cigar}.')
            if tx_id not in transcripts:
                transcripts[tx_id] = {
                    'gen_chrom': chrom,
                    'gen_start': int(pos),
                    'cigar': cigar
                }
            else:
                raise ValueError(f'Transcript {tx_id} may be duplicated.')
    return transcripts


def is_cigar_valid(cigar_str: str) -> bool:
    '''
    Checks for a valid CIGAR string. 
    "Valid" means containing a number 
    followed by one of X (mismatch), I (insertion), D (deletion), or M (match). CIGAR
    strings containing S, H, or = are not covered by this program. 
    Regex uses * for zero or more times, as 
    "in some CIGAR variants, the integer may be omitted if it is 1."
    (from https://www.drive5.com/usearch/manual/cigar.html)  
    
    Args:
        cigar_str: cigar string to validate
    Returns:
        Bool (True or False) based on CIGAR validity 
    '''

    return bool(re.match(r'^([0-9]*[M|I|D|X])+$', cigar_str))


def convert_cigar_string_to_list(cigar_str: str) -> list:
    '''
    Uses RegEx to find all instances of a number (\\d) repeated >=1 times followed by 
    a word character and return that list. 

    Uses RegEx to look to see if there are two letters consecutively, in which case, 
    re.sub will add a 1 between them (e.g., 7DM --> 7D1M). re.sub is not currently
    configured to handle the case in which there are >2 letters in a row (like IDXM).

    Args: 
        cigar_str: cigar string to convert

    Returns:
        list of tuples (span, cigar_type):
        8M7D6X2I2M11D7M --> 
        [('8', 'M'),
         ('7', 'D'),
         ('6', 'X'),
         ('2', 'I'),
         ('2', 'M'),
         ('11', 'D'),
         ('7', 'M')]
    '''

    mod_cig = re.sub(r'([DIMX])([DIMX])', r'\g<1>1\2', cigar_str)
    cigar_list = re.findall(r'(\d+)(\w)', mod_cig)
    return cigar_list


def get_genome_pos(tr_coord: int, cigar: str, start_pos: int) -> int:
    '''
    Takes transcript coordinate, CIGAR string, and genome start position 
    and calculates corresponding genome_position of transcript coordinate.

    Args:
        tr_coord: the coordinate of interest/of the query
        cigar: valid CIGAR string (how the transcript maps to the genome)
        start_pos: genome start position

    Returns:
        gen_pos: genomic position that corresponds to transcript coordinate
    '''
    gen_pos = int(start_pos) - 1 #adjusting for 0-based index
    tr_remain = int(tr_coord) + 1 #adjusting for 0-based index
    cigar_list = convert_cigar_string_to_list(cigar)

    for span, cigar_type in cigar_list:
        #if there are no more coords remaining, we're done!
        if int(tr_remain) == 0:
            break
        if cigar_type in ('M', 'X'):
            #if it's a match or mismatch, both the gen_pos and tr_remain numbers are affected
            #(e.g., we advance through both sequences)
            #GENOME:CHR1 ACTGTCATGTACGTTTAGCTAGCC--TAGCTAGGGACCTAGATAATTTAGCTAG
            #TR1            GTCATGTA-------CTAGCCGGTA-----------AGATAAT
            gen_pos += min(int(span), tr_remain)
            tr_remain -= min(tr_remain, int(span))
        elif cigar_type == 'D':
            #if there is a deletion in the transcript, only the gen_pos is affected
            gen_pos += int(span)
        elif cigar_type == 'I':
            #if there is an insertion in the transcript,
            #only the tr_remain is affected (see above drawing)
            tr_remain -= min(tr_remain, int(span))
    if tr_remain > 0:
        #after this for loop, if there are still bases remaining, then it is not a
        #valid query.
        print(f'Found invalid query starting from genomic coord {gen_pos}, \
            using {cigar} CIGAR -- will be returned as -1 in output file!')
        return -1
    return gen_pos

def query_transcript(transcript_file, query_file, output_file) -> None:
    '''
    Writes output file such that for every query in query_file, each line contains 4 cols:
    query transcript, query transcript position, the genomic
    chromosome, and the genomic position of the query position.

    Args:
        Transcript_file: path to the transcript file (.tsv) -- each line is one transcript
            composed of 4 columns: transcript ID, genomic chromosome, 
            genomic position, and CIGAR string.
        Query_file: path to the query file (.tsv) -- each line is one query composed of
            2 columns: transcript ID of query, and transcript coordinate of query
        Output_file: path to output tsv -- for each query in query_file, 
            each line in output_file composed of 4 columns: transcript ID,
            transcript coordinate, genomic chromosome, and genomic position.
    Returns:
        --
    Raises:
        ValueError: if either of the input files are malformed
    '''

    # Build transcript dict
    transcripts = assemble_transcript_dict(transcript_file=transcript_file)

    # Iterate through queries
    with open(output_file, 'w', encoding='utf-8') as out_f:
        out_writer = csv.writer(out_f, delimiter='\t') #make sure output is tab-separated
        with open(query_file, 'r', encoding='utf-8') as query_f:
            for line in query_f:
                fields = [i.strip() for i in line.split('\t')]
                if len(fields) != 2:
                    raise ValueError(f'Unexpected format of query file, \
                        expected 2 columns, got {len(fields)} columns.')
                tx_id, tr_coord = fields[0], int(fields[1])
                if tr_coord < 0:
                    raise ValueError(f'Coord {tr_coord} in {tx_id} is negative, and thus is not considered valid input!')
                # Get query transcript info from transcript dict
                gen_chrom = transcripts[tx_id]['gen_chrom']
                gen_start = transcripts[tx_id]['gen_start']
                cigar_str = transcripts[tx_id]['cigar']
                gen_pos = get_genome_pos(
                    tr_coord=tr_coord,
                    cigar=cigar_str,
                    start_pos=gen_start)
                out_writer.writerow([tx_id, tr_coord, gen_chrom, gen_pos])

def parse_args():
    '''
    Parse arguments during runtime.
    Args:
        --
    Returns:
        User-specified arguments for transcript_fie, query_file, and out.
    Raises:
        --
    '''
    parser = ArgumentParser(description='This script translates transcript \
                            coordinates to genomic coordinates')
    parser.add_argument('-t', '--transcript-file', type=str, required=True,
                        help='Transcript file, 4 columns required:, name, CHR, POS, CIGAR')
    parser.add_argument('-q', '--query-file', type=str, required=True,
        help='Query file, 2 columns required: transcript name, transcript coord')
    parser.add_argument('-o', '--out', type=str, required=True,
        help='Output file to write results to')

    args = parser.parse_args()
    return args.transcript_file, args.query_file, args.out


def main(transcript_file, query_file, output_file):
    '''
    Args:
        transcript_file: path to the transcript file (.tsv)
        query_file: path to query file (.tsv)
        out: path to output file (.tsv).
  
    Returns:
        --
    '''
    query_transcript(
        transcript_file=transcript_file,
        query_file=query_file,
        output_file=output_file)

if __name__ == '__main__':
    main(
        *parse_args()
    )
