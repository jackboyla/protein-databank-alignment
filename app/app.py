
#!/bin/python3

import math
import os
import random
import re
import sys

#
# Complete the 'find_pdb_mapping' function below.

# The function accepts following arguments:
#  1. (str) reference_sequence: the primary structure of the protein
#  2. (List[Tuple[str, str]]) pdb_sequence: list of (PDB position, amino acid 1 letter code) tuples.

# PDB positions are represented as strings with the format "{chain_id}.{position}.{insertion_code}"
# For example:
# pos_1 = "A.1."
# pos_2 = "A.2."
# pos_3 = "A.2.A"

# The function is expected to return a list of Tuple[str, int] linking pdb positions to indexes (starting from 0)
# in the reference sequence.

# For instance: [("A.1.", 0), ("A.2.", 1), ("A.2.A", 2)]
#

def find_pdb_mapping(reference_sequence, pdb_sequence):
    '''
    inputs:
        reference_sequence [str] e.g, "MMAAGHSFDJNDKVDKNKNDLNDFLVKNDFKNVKNXFVSSDNDOSRJIODRJ"
        pdb_sequence [List] e.g, [('A.254.', 'A'), ('A.255.', 'G'), ('A.258.', 'F'), ('A.258.A', 'D'), ('A.264.', 'D')]

    returns: [('A.254.', 2), ('A.255.', 4), ('A.258.', 7), ('A.258.A', 8), ('A.264.', 14)]

    NOTE:   This implementation assumes there is only ONE correct match in the reference sequence.

            Also, this assumes there may be N insertion aminos between consecutive locations
            i.e, between A.254 and A.255
            If this is not the case, then the regex building must be changed to reflect this.

            Also, this assumes that insertion codes must be consecutively listed ('', 'A', 'B', 'C', ...) and
            that all previous insertion codes for that location must be listed, i.e,
                'A.258.A' cannot be listed unless 'A.258.' is also listed.
            The only exception is if the location is the first entry in the pdb sequence.
    '''
    
    # combine insertion aminos and update idx to reflect
    # the minimum gap we can expect between aminos
    # eg, [('A.254.', 'A'), ('A.255.', 'G'), ('A.258.', 'F'), ('A.258.A', 'D'), ('A.264.', 'D')]
    # ->  [(254, 'A'), (255, 'G'), (258, 'F', 'D'), (265, 'D')]
    seq_idx = parse_pdb_seq(pdb_sequence)

    reg, gaps = get_regex_with_gaps(seq_idx)


    # reg = A[A-Z]{0,}?G[A-Z]{2,}?FD[A-Z]{5,}?D

    # find the section of the ref sequence that satisfies the pdb sequence
    res = re.search(reg, reference_sequence)  # AGHSFDJNDKVD
    match_start = res.start()
    match_idx = [(i + match_start, c) for i, c in enumerate(res.group())]

    out = align_match_to_pdb_seq(match_idx, pdb_sequence, gaps)

    return out


def convert_loc_to_idx(pdb_tuple):
    pos = pdb_tuple[0].split(".")
    idx = int(pos[1])
    return idx


def parse_pdb_seq(pdb_sequence):
    seq_idx = []
    for pos, amino in pdb_sequence:
        idx = convert_loc_to_idx((pos, amino))
        if pos.split(".")[-1] != "":
            if len(seq_idx) > 0:
                seq_idx[-1] += (amino,)
            else:
                seq_idx.append((idx, amino))
        else:
            seq_idx.append((idx, amino))
    
    return seq_idx


def get_regex_with_gaps(seq_idx):
    # build regex for finding the pdb_sequence in the ref seq
    # and observe the minimum gaps between aminos
    reg = r""
    gaps = []
    for i, tup in enumerate(seq_idx):

        if i == 0:
            for amino in tup[1:]:
                reg += re.escape(amino)
                gaps.append(0)

        else: 
            # we don't know exactly how many aminos may be between position 1 and 2,
            # (because there could be N insertion aminos)
            # but we know it must be at least 1
            gap_min = tup[0] - seq_idx[i-1][0] -1 
            gaps.append(gap_min)

            # if the aminos are in an insertion group, 
            # they come directly after each other
            if len(tup) > 2:
                exact = re.escape(tup[1])
                for amino in tup[2:]:
                    exact += re.escape(amino)
                    gaps.append(0)
            else:
                exact = re.escape(tup[1])
        
            # we add the "?" due to *greediness*: https://stackoverflow.com/a/16619508
            # we want the regex to search from shortest to longest match, not longest to shortest
            reg += r"[A-Z]{" + re.escape(str(gap_min)) + r",}?" + re.escape(exact)

    return reg, gaps


def align_match_to_pdb_seq(match_idx, pdb_sequence, gaps):

    # iterate through the matched string to find exact locations of the pdb_sequence aminos
    k = 0
    out = [(pdb_sequence[k][0], match_idx[k][0])]

    i = k = 1
    while len(out) != len(pdb_sequence):

        # if the aminos match...
        if match_idx[i][1] == pdb_sequence[k][1]:

            # if the gap is sufficiently large...
            if (match_idx[i][0] - out[k-1][1] ) > (gaps[k]):
                out.append((pdb_sequence[k][0], match_idx[i][0]))
                k += 1
        i+=1
    return out
    
    
            
    

if __name__ == '__main__':

    reference_sequence = "MMAAGHSFDJNDKVDKNKNDLNDFLVKNDFKNVKNXFVSSDNDOSRJIODRJ"
    pdb_sequence = [('A.254.', 'A'), ('A.255.', 'G'), ('A.258.', 'F'), ('A.258.A', 'D'), ('A.264.', 'D')]

    out = find_pdb_mapping(reference_sequence, pdb_sequence)
    print(out)
    # fptr = open(os.environ['OUTPUT_PATH'], 'w')

    # reference_sequence = input()

    # pdb_sequence_count = int(input().strip())

    # pdb_sequence = []

    # for _ in range(pdb_sequence_count):
    #     pdb_sequence_item = input()
    #     position, amino_acid = pdb_sequence_item.split(" ")
    #     pdb_sequence.append((position, amino_acid))

    # pdb_mapping_result = find_pdb_mapping(reference_sequence, pdb_sequence)

    # result = [position + " " + str(index) for position, index in pdb_mapping_result]
    # fptr.write('\n'.join(result))
    # fptr.write('\n')

    # fptr.close()
