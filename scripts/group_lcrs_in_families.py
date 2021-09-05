import numpy as np
import os
import align_mafft as am
import get_lcrs_stats as stats

from typing import Dict, List, Tuple


def get_lcrs_from_maffted(mafft_output: Dict[stats.Header, str], lcrs: Dict[stats.PfamHeaderLCR, str]) -> List[Tuple[str, int, int, str, str]]:
    lcrs_from_maffted = []
    for header in mafft_output:
        seq = mafft_output[header]
        for lcr in lcrs:
            if lcr.get_entry_name() == header.get_entry_name():
                start = lcr.get_lcr_start()
                end = lcr.get_lcr_end()
                num_dashes = 0
                is_lcr = False
                maffted_lcr = ""
                for i, aa in enumerate(seq):
                    if i == start and aa != "-":
                        is_lcr = not is_lcr
                    elif i == end and aa != "-" and is_lcr:
                        maffted_lcr += aa
                        lcrs_from_maffted.append((lcrs[lcr], start, end + num_dashes, lcr.__str__(), maffted_lcr))
                        break
                    if is_lcr:
                        maffted_lcr += aa
                    if aa == "-":
                        start += 1
                        end += 1
                        if is_lcr:
                            num_dashes += 1
    print(sorted(lcrs_from_maffted, key=lambda e: len(e[0]), reverse=True))
    return sorted(lcrs_from_maffted, key=lambda e: len(e[0]), reverse=True)


def count_not_common(s1: int, e1: int, m1: str, s2: int, e2: int, m2: str) -> int:
    if s1 > s2:
        print('p1')
        m2 = m2[s1-s2:]
    elif s1 < s2:
        print('p2')

        m1 = m1[s2-s1:]
    if e1 > e2:
        print('p3')

        m1 = m1[:e2-e1]
    elif e1 < e2:
        print('p4')
        print(e1-e2)
        print(m1)
        print(m2)
        m2 = m2[:e1-e2]
        print(m2)

    if len(m1) != len(m2):
        print(s1, e1, m1)
        print(s2, e2, m2)
        raise Exception('Method is wrong - common parts have not same length!')

    dashes = 0
    for i in range(len(m1)):
        if m1[i] != '-' and m2[i] != '-':
            dashes += 1
    return dashes


def write_grouped_lcrs(lcrs_from_mafft: List[Tuple[str, int, int, str, str]], output_name: str, family_name: str):
    file_name = "./{}/{}.fasta".format(output_name, family_name)
    if os.path.dirname(file_name):
        os.makedirs(os.path.dirname(temp_file_name), exist_ok=True)
    file = open(file_name, "w")
    while lcrs_from_mafft:
        file.write(">" + lcrs_from_mafft[0][3] + " | maffted: " + str(lcrs_from_mafft[0][1]) + " - " + str(lcrs_from_mafft[0][2]) + "\n")
        file.write(lcrs_from_mafft[0][0] + "\n")

        for lcr in lcrs_from_mafft[1:]:
            if count_not_common(lcrs_from_mafft[0][1], lcrs_from_mafft[0][2], lcrs_from_mafft[0][4], lcr[1], lcr[2], lcr[4]) >= parameter * len(lcr[0]):
                file.write(">" + lcr[3] + " | maffted: " + str(lcr[1]) + " - " + str(lcr[2]) + "\n")
                file.write(lcr[0] + "\n")
                lcrs_from_mafft.remove(lcr)
        lcrs_from_mafft = lcrs_from_mafft[1:]
        file.write("-----\n")
    file.close()


if __name__ == "__main__":
    # lcrs_file_name = "./../LCR_sequences/Pfam_LCRs.fasta"
    lcrs_file_name = "./../LCR_sequences/test_mafft_lcrs.fasta"
    header_type = stats.PfamHeaderLCR
    # pfam_file_name = "./../sequences/Pfam-A.fasta"
    pfam_file_name = "./../LCR_sequences/test_mafft.fasta"
    # output_name = "mafft_Pfam_LCRs_grouped"
    output_name = "mafft_Pfam_LCRs_test2"
    parameter = 0.5
    lcr_families = am.get_full_seqs_per_family(stats.get_seqs_dict(lcrs_file_name, header_type))
    mafft_lcr_families = {k: v for k, v in lcr_families.items()
                          if len(v) > 1 and len(set([i.get_identification_str() for i in v])) > 1}
    prot_seqs = am.get_common_pfam_seqs_per_family(pfam_file_name, mafft_lcr_families)

    for family in mafft_lcr_families:
        family_prot_seqs = {header: prot_seqs[header] for header in prot_seqs
                            if header.__str__() in [i.__str__().split(" LCR")[0]
                                                    for i in mafft_lcr_families[family]]}
        temp_file_name = am.write_temp_file(family_prot_seqs, output_name)
        mafft_file_name = am.run_mafft(temp_file_name, family, output_name)
        os.remove(temp_file_name)

        lcrs_from_maffted = get_lcrs_from_maffted(am.get_mafft_output(mafft_file_name, stats.PfamHeader),
                                                  mafft_lcr_families[family])
        write_grouped_lcrs(lcrs_from_maffted, output_name, family)
        os.remove(mafft_file_name)
