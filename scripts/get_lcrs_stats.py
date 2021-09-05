"""
Script by Maciej Dzikowski.
"""

import numpy as np
import matplotlib.pyplot as plt

from abc import ABC, abstractmethod
from typing import Callable, Dict


class Header:
    def __init__(self, raw_header: str):
        self.header = raw_header.strip().strip(">")

    def __str__(self):
        return self.header

    def get_one_word_info(self, abbr: str) -> str:
        return self.header.split(abbr + "=")[1].split()[0]


class UniprotHeader(Header):
    def get_database_id(self) -> str:
        return self.header.split("|")[0]

    def get_unique_identifer(self) -> str:  # uniprot_id
        return self.header.split("|")[1]

    def get_entry_name(self) -> str:
        return self.header.split()[0].split("|")[2]

    def get_protein_name(self) -> str:
        return " ".join(self.header.split(" OS=")[0].split()[1:])

    def get_organism_name(self) -> str:
        return self.header.split("OS=")[1].split(" OX=")[0]

    def get_organism_identifier(self) -> str:
        return self.get_one_word_info("OX")

    def get_gene_name(self) -> str:
        if "GN=" in self.header:
            return self.get_one_word_info("GN")

    def get_protein_existence(self) -> int:
        return int(self.get_one_word_info("PE"))

    def get_sequence_version(self) -> int:
        return int(self.get_one_word_info("SV"))


class PfamHeader(Header):
    def get_entry_name(self) -> str:
        return self.header.split("/")[0]

    def get_prot_start(self) -> int:
        return int(self.header.split("/")[1].split("-")[0])

    def get_prot_end(self) -> int:
        return int(self.header.split("-")[1].split()[0])

    def get_unique_identifer(self) -> str:  # uniprot_id + "." + <int?>
        return self.header.split()[1]

    def get_family_id(self) -> str:
        return self.header.split()[2].split(".")[0]

    def get_family_name(self) -> str:
        return self.header.split()[2].split(";")[1]

    def get_identification_str(self) -> str:
        return self.get_entry_name() + "/" + str(self.get_prot_start()) + "-" + str(self.get_prot_end())

class LCRHeader(Header, ABC):
    def get_lcr_start(self) -> int:
        return int(self.get_one_word_info("LCR:begin").strip(","))

    def get_lcr_end(self) -> int:
        return int(self.get_one_word_info("end"))


class UniprotHeaderLCR(UniprotHeader, LCRHeader):
    pass


class PfamHeaderLCR(PfamHeader, LCRHeader):
    pass


class PfamInSPHeaderLCR(PfamHeaderLCR, UniprotHeaderLCR):
    pass


class Statistics(ABC):
    def __init__(self, file_name: str):
        self.file_name = file_name

        self.seqs = get_seqs_dict(file_name, self.get_header())

    @abstractmethod
    def get_header(self) -> Callable:
        pass

    def get_num_seqs(self) -> int:
        return len(self.seqs)

    def get_seqs_lengths(self) -> list:
        return [len(seq) for seq in self.seqs.values()]

    def get_min_max_avg_sd(self, data: list) -> tuple:
        return min(data), max(data), np.mean(data), np.std(data)

    @abstractmethod
    def print_statistics(self):
        pass


class UniprotStatistics(Statistics):
    def get_header(self) -> Callable:
        return UniprotHeader

    def get_seqs_per_organism(self) -> dict:
        organisms = {}
        for seq in self.seqs:
            org = seq.get_organism_name()
            if org in organisms:
                organisms[org] += [self.seqs[seq]]
            else:
                organisms[org] = [self.seqs[seq]]
        return organisms

    def get_num_organisms(self) -> int:
        return len(self.get_seqs_per_organism())

    def get_num_seqs_per_organism(self) -> list:
        organisms = self.get_seqs_per_organism()
        return [len(organisms[org]) for org in organisms]

    def get_avg_len_per_organism(self) -> list:
        organisms = self.get_seqs_per_organism()
        return [np.mean([len(seq) for seq in organisms[org]]) for org in organisms]

    def print_statistics(self):
        print("Statistics for", self.file_name)
        print("Number of sequences:", self.get_num_seqs())
        print("Length of sequences:", self.get_min_max_avg_sd(self.get_seqs_lengths()))
        print("Number of organisms:", self.get_num_organisms())
        print("Number of sequences per organism:", self.get_min_max_avg_sd(self.get_num_seqs_per_organism()))
        print("Length of sequences per organism:", self.get_min_max_avg_sd(self.get_avg_len_per_organism()))
        print("---", "\n")


class PfamStatistics(Statistics):
    def get_header(self) -> Callable:
        return PfamHeader

    def get_seqs_per_family(self) -> dict:
        families = {}
        for seq in self.seqs:
            family = seq.get_family_name()
            if family in families:
                families[family] += [self.seqs[seq]]
            else:
                families[family] = [self.seqs[seq]]
        return families

    def get_num_families(self) -> int:
        return len(self.get_seqs_per_family())

    def get_num_seqs_per_family(self) -> list:
        families = self.get_seqs_per_family()
        return [len(families[family]) for family in families]

    def get_avg_len_per_family(self) -> list:
        families = self.get_seqs_per_family()
        return [np.mean([len(seq) for seq in families[family]]) for family in families]

    def print_statistics(self):
        print("Statistics for", self.file_name)
        print("Number of sequences:", self.get_num_seqs())
        print("Length of sequences:", self.get_min_max_avg_sd(self.get_seqs_lengths()))
        print("Number of families:", self.get_num_families())
        print("Number of sequences per family:", self.get_min_max_avg_sd(self.get_num_seqs_per_family()))
        print("Length of sequences per family:", self.get_min_max_avg_sd(self.get_avg_len_per_family()))
        print("---", "\n")


class UniprotStatisticsLCR(UniprotStatistics):
    def get_header(self) -> Callable:
        return UniprotHeaderLCR

    def get_seqs_per_prot(self) -> dict:
        prots = {}
        for seq in self.seqs:
            prot = seq.get_entry_name()
            if prot in prots:
                prots[prot] += [self.seqs[seq]]
            else:
                prots[prot] = [self.seqs[seq]]
        return prots

    def get_num_prots(self) -> int:
        return len(self.get_seqs_per_prot())

    def get_num_seqs_per_prot(self) -> list:
        prots = self.get_seqs_per_prot()
        return [len(prots[prot]) for prot in prots]

    def get_avg_len_per_prot(self) -> list:
        prots = self.get_seqs_per_prot()
        return [np.mean([len(seq) for seq in prots[prot]]) for prot in prots]

    def get_seqs_per_org(self) -> dict:
        orgs = {}
        for seq in self.seqs:
            org = seq.get_protein_name()
            if org in orgs:
                orgs[org] += [self.seqs[seq]]
            else:
                orgs[org] = [self.seqs[seq]]
        return orgs

    def get_num_orgs(self) -> int:
        return len(self.get_seqs_per_org())

    def get_num_seqs_per_org(self) -> list:
        orgs = self.get_seqs_per_org()
        return [len(orgs[org]) for org in orgs]

    def get_avg_len_per_org(self) -> list:
        orgs = self.get_seqs_per_org()
        return [np.mean([len(seq) for seq in orgs[org]]) for org in orgs]

    def print_statistics(self):
        print("Statistics for", self.file_name)
        print("Number of sequences:", self.get_num_seqs())
        print("Length of sequences:", self.get_min_max_avg_sd(self.get_seqs_lengths()))
        print("Number of prot. sequences with LCRs:", self.get_num_prots())
        print("Length of prot. sequences with LCRs:", )
        print("Number of LCRs per prot. sequence:", )
        print("Number of LCRs per prot. sequence with LCR:", self.get_min_max_avg_sd(self.get_num_seqs_per_prot()))
        print("Length of LCRs per prot. sequence:")
        print("Length of LCRs per prot. sequence with LCR:", self.get_min_max_avg_sd(self.get_avg_len_per_prot()))
        print("Number of organisms:", self.get_num_orgs())
        print("Number of organism with LCRs:", )
        print("Length of organism with LCRs:", )
        print("Number of LCRs per organism:", )
        print("Number of LCRs per organism with LCR:", self.get_min_max_avg_sd(self.get_num_seqs_per_org()))
        print("Length of LCRs per organism:", )
        print("Length of LCRs per organism with LCR:", self.get_min_max_avg_sd(self.get_avg_len_per_org()))
        print("---", "\n")


class PfamStatisticsLCR(PfamStatistics):
    def get_header(self) -> Callable:
        return PfamHeaderLCR

    def get_seqs_per_prot(self) -> dict:
        prots = {}
        for seq in self.seqs:
            prot = seq.get_entry_name()
            if prot in prots:
                prots[prot] += [self.seqs[seq]]
            else:
                prots[prot] = [self.seqs[seq]]
        return prots

    def get_num_prots(self) -> int:
        return len(self.get_seqs_per_prot())

    def get_num_seqs_per_prot(self) -> list:
        prots = self.get_seqs_per_prot()
        return [len(prots[prot]) for prot in prots]

    def get_avg_len_per_prot(self) -> list:
        prots = self.get_seqs_per_prot()
        return [np.mean([len(seq) for seq in prots[prot]]) for prot in prots]

    def print_statistics(self):
        print("Statistics for", self.file_name)
        print("Number of LCRs:", self.get_num_seqs())
        print("Length of LCRs:", self.get_min_max_avg_sd(self.get_seqs_lengths()))
        print("Number of prot. sequences with LCRs:", self.get_num_prots())
        print("Length of prot. sequences with LCRs:", )
        print("Number of LCRs per prot. sequence:", )
        print("Number of LCRs per prot. sequence with LCR:", self.get_min_max_avg_sd(self.get_num_seqs_per_prot()))
        print("Length of LCRs per prot. sequence:")
        print("Length of LCRs per prot. sequence with LCR:", self.get_min_max_avg_sd(self.get_avg_len_per_prot()))
        print("Number of families:", self.get_num_families())
        print("Number of family with LCRs:", )
        print("Length of family with LCRs:", )
        print("Number of LCRs per family:", )
        print("Number of LCRs per family with LCR:", self.get_min_max_avg_sd(self.get_num_seqs_per_family()))
        print("Length of LCRs per family:", )
        print("Length of LCRs per family with LCR:", self.get_min_max_avg_sd(self.get_avg_len_per_family()))
        print("---", "\n")


class PfamInSPStatisticsLCR(UniprotStatistics, PfamStatistics):
    def get_header(self) -> Callable:
        return PfamInSPHeaderLCR

    def print_statistics(self):
        print("Statistics for", self.file_name)
        print("Number of sequences:", self.get_num_seqs())
        print("Length of sequences:", self.get_min_max_avg_sd(self.get_seqs_lengths()))
        print("---", "\n")


def get_seqs_dict(file_name: str, header_type: Callable) -> Dict[Header, str]:
    seqs = {}
    header = False
    with open(file_name) as f:
        for line in f:
            if line.startswith(">"):
                header = header_type(line)
                seqs[header] = ""
            elif line:
                if header:
                    seqs[header] += line.strip()
    return seqs


def check_dbs_identical(db1: Statistics, db2: Statistics) -> bool:
    checker = True
    if len(db1.seqs) == len(db2.seqs):
        db2_seqs = {seq.header: "" for seq in db2.seqs}
        for seq in db1.seqs:
            if seq.header not in db2_seqs.keys():
                checker = False
                print(seq.header)
    else:
        checker = False
    return checker


if __name__ == "__main__":
    # UniprotStatistics('./../sequences/uniprot_sprot.fasta').print_statistics()
    # PfamStatistics('./../sequences/Pfam-A.fasta').print_statistics()
    # UniprotStatisticsLCR('./../LCR_sequences/swiss_prot_LCRs_052021.fasta').print_statistics()
    # PfamStatisticsLCR('./../LCR_sequences/Pfam_LCRs_intermediate.fasta').print_statistics()
    #
    # pfam_in_sp = PfamInSPStatisticsLCR('./../LCR_sequences/Pfam-sp_052021-LCRs.fasta')
    # pfam_in_sp_lcrs = PfamInSPStatisticsLCR('./../LCR_sequences/Pfam-sp_LCRs_052021-LCRs.fasta')
    # if check_dbs_identical(pfam_in_sp, pfam_in_sp_lcrs):
    #     pfam_in_sp.print_statistics()
    # else:
    #     pfam_in_sp.print_statistics()
    #     pfam_in_sp_lcrs.print_statistics()

    lens = UniprotStatisticsLCR('./../LCR_sequences/swiss_prot_LCRs_052021.fasta').get_seqs_lengths()
    x = [i for i in range(1, 2034)]
    y = [lens.count(i+1) for i in range(2033)]

    # plt.plot(x, y, 'k')
    # plt.xlabel("length in amino acids")
    # plt.ylabel("number of LCRs")
    # plt.savefig("./len_dist.png")

    # plt.plot(x, y, 'k')
    # plt.xlim(1, 100)
    # plt.savefig("./len_dist_100.png")

    # fig, (ax1, ax2) = plt.subplot(2)
    # plt.ylim(0, 10)
    # plt.savefig("./len_dist_y.png")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.set_title('$x \in <1, 100>$')
    ax1.set_xlabel("length in amino acids")
    ax1.set_ylabel("number of LCRs")
    ax1.set_xlim(1, 100)
    ax1.plot(x, y, color='black')

    ax2.set_title('$y \in <0, 10>$')
    ax2.set_xlabel("length in amino acids")
    ax2.set_ylim(0, 10)
    ax2.plot(x, y, color='black')

    plt.subplots_adjust(hspace=1.5)
    plt.savefig("./len_dist_sub.png")
