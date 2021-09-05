import re
import get_lcrs_stats as stats

from typing import List, Set


class Sequence:
    def __init__(self, header: stats.Header, seq: str):
        self.header = header
        self.seq = seq

    def __str__(self):
        return self.header.__str__()


def create_seq_dicts(file: str) -> List[Sequence]:
    seq_list = []
    header = False
    with open(file) as f1:
        for line in f1:
            if line.startswith(">"):
                if header:
                    seq_list.append(Sequence(header, seq))
                header = stats.Header(line.strip(';\n'))
                seq = ""
            elif line:
                if header:
                    seq += line.strip()
    return seq_list


def get_matrix_keys(matrix_file: str) -> Set[str]:
    matrix_dict = {}
    amino_acids_list = []
    for line in open(matrix_file):
        if not line.startswith("#") and re.match(r'[ \t]', line):
            amino_acids_list = line.split()
    print(amino_acids_list)
    return set(amino_acids_list)


if __name__ == "__main__":
    sequences = create_seq_dicts("./../LCR_sequences/swiss_prot_LCRs_052021.fasta")
    keys = get_matrix_keys("./../prob_matrices/exp_prob_matrix")

    print("Szukanie...")
    to_trash = []
    for seq in sequences:
        for aa in seq.seq:
            if aa not in keys:
                to_trash.append(seq)
                break

    for i in to_trash:
        print(i.header, i.seq)
