Goal: Translate a (0-based) transcript coordinate to a (0 based) genome coordinate 
given the following information: chromosome, genomic start position, CIGAR string 
(which indicates how the transcript maps to the genome sequence), and the transcript
coordinate to query.

1. Assumptions:
    * Transcript is always mapped from genomic 5' -> 3'
    * Do not need to accomodate rarer CIGAR annotations (S, H, =, etc.)
    * Empty CIGAR strings are not valid input 
    * Tab-separated input is provided  

2. Strengths:
    * Checks for valid input (number of columns in input files, CIGAR string formation)
    * Run time dependent on CIGAR length but function stops when tr_coord of interest is reached
    (so we don't have to go through the whole CIGAR string) 

3. Weaknesses:

    * My solution does not accomodate the edge case where the CIGAR has >2 letters in a row 
    * Line-by-line parsing does not scale well; could use pandas.read_csv, 
    which is much faster than csvreader under the hood. 
    * RegEx adds some additional overhead; could just iterate through the CIGAR string as needed
    to ensure valid formatting


4. Improvements
    * Speed: use pandas to read-in data, refactor to use .applymap() with get_genome_pos()
