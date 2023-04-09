import unittest

from app.app import find_pdb_mapping


class PDBMappingTest(unittest.TestCase):

    def test_case(self):
        expected_result = [
            [('A.254.', 2), ('A.255.', 4), ('A.258.', 7), ('A.258.A', 8), ('A.264.', 14)],
            [('B.267.', 5), ('B.267.A', 6), ('B.267.B', 7), ('B.267.C', 8), ('B.270.', 11)],
            [('C.140.B', 0), ('C.140.C', 1), ('C.140.D', 2), ('C.143.', 6), ('C.145.', 9)],
            [('D.34.C', 1), ('D.34.D', 2), ('D.34.D', 3), ('D.34.E', 4), ('D.35.', 5), ('D.35.A', 6), ('D.36.', 9), ('D.40.', 27)]
                           ]

        reference_sequence = [
            "MMAAGHSFDJNDKVDKNKNDLNDFLVKNDFKNVKNXFVSSDNDOSRJIODRJ",
            "SDAAAAAFNVLKDKVMKNVLSDNLVJSDNVKJSNVKLJKNDVLKSNDLVNSIGUI",
            "KSDNBJNOPODPDPOSNCKLMGPNF",
            "LKSDNVOEBOICSKJDNOJEVNSDKLCX"
                              ]
        
        pdb_sequence = [
            [('A.254.', 'A'), ('A.255.', 'G'), ('A.258.', 'F'), ('A.258.A', 'D'), ('A.264.', 'D')],
            [('B.267.', 'A'), ('B.267.A', 'A'), ('B.267.B', 'F'), ('B.267.C', 'N'), ('B.270.', 'K')],
            [('C.140.B', 'K'), ('C.140.C', 'S'), ('C.140.D', 'D'), ('C.143.', 'N'), ('C.145.', 'O')],
            [('D.34.C', 'K'), ('D.34.D', 'S'), ('D.34.D', 'D'), ('D.34.E', 'N'), ('D.35.', 'V'), ('D.35.A', 'O'), ('D.36.', 'O'), ('D.40.', 'X')]
                        ]


        for r, p, e in zip(reference_sequence, pdb_sequence, expected_result):

            result = find_pdb_mapping(r, p)
    
            self.assertEqual(result, e)
