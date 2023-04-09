Complete the 'find_pdb_mapping' function in app.py.

The function accepts following arguments:
 1. (str) reference_sequence: the primary structure of the protein
 2. (List[Tuple[str, str]]) pdb_sequence: list of (PDB position, amino acid 1 letter code) tuples.

PDB positions are represented as strings with the format "{chain_id}.{position}.{insertion_code}"
For example:
pos_1 = "A.1."
pos_2 = "A.2."
pos_3 = "A.2.A"

The function is expected to return a list of Tuple[str, int] linking pdb positions to indexes (starting from 0)
in the reference sequence.

For instance: [("A.1.", 0), ("A.2.", 1), ("A.2.A", 2)]


    NOTE:   This implementation assumes there is only ONE correct match in the reference sequence.

            Also, this assumes there may be N insertion aminos between consecutive locations
            i.e, between A.254 and A.255
            If this is not the case, then the regex building must be changed to reflect this.

            Also, this assumes that insertion codes must be consecutively listed ('', 'A', 'B', 'C', ...) and
            that all previous insertion codes for that position must be listed, i.e,
                'A.258.A' cannot be listed unless 'A.258.' is also listed