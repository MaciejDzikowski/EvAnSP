import re
import scripts.get_lcrs_stats as stats

from typing import Dict, List, Tuple


class Sequence:
    def __init__(self, header: stats.Header, seq: str):
        self.header = header
        self.seq = seq

    def __str__(self):
        return self.header.__str__()


class Parser:
    def __init__(self, first_seq: str, second_seq: str, matrix_file: str):
        self.first_seq = first_seq
        self.second_seq = second_seq
        self.matrix_file = matrix_file

        self.first_seq_dict = create_seq_list(self.first_seq)
        self.second_seq_dict = create_seq_list(self.second_seq)
        self.matrix = self.create_matrix()

    def create_matrix(self) -> Dict[str, float]:
        matrix_dict = {}
        amino_acids_list = []
        for line in open(self.matrix_file):
            if not line.startswith("#") and re.match(r'[ \t]', line):
                amino_acids_list = line.split()
            if not line.startswith("#") and not re.match(r'[ \t]', line):
                if not amino_acids_list:
                    raise Exception('no amino_acids_list')
                line_list = line.split()
                for i, num in enumerate(line_list[1:]):
                    matrix_dict[line_list[0] + amino_acids_list[i]] = float(num)
        return matrix_dict

    def get_parsered_input(self) -> Tuple[List[Sequence], List[Sequence], Dict[str, float]]:
        print("Parsing data...")
        return self.first_seq_dict, self.second_seq_dict, self.matrix

    def check_file(self):
        return None


def create_seq_list(file: str) -> List[Sequence]:
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


def get_evansp_output(file: str) -> Tuple[list, list]:
    results_list = []
    same_results_list = []
    counter = 0
    with open(file) as f:
        seq1 = None
        seq2 = None
        prob = None
        for line in f:
            if counter % 5 == 0:
                seq1 = line.strip()
            if counter % 5 == 1:
                seq2 = line.strip()
            if counter % 5 == 3:
                prob = float(line.strip().split(" frame")[0].split()[1])
                results_list.append(prob)
                if seq1 == seq2:
                    same_results_list.append(prob)
            counter += 1
    return results_list, same_results_list


def sort_evansp_output(file: str) -> dict:
    pass


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    file_names = ["./../results/EvAnSP_01092021125511/2_EvAnSP_01092021125511.fasta",
                  "./../results/EvAnSP_01092021130312/3_EvAnSP_01092021130312.fasta",
                  "./../results/EvAnSP_01092021131251/4_EvAnSP_01092021131251.fasta",
                  "./../results/EvAnSP_01092021132420/5_EvAnSP_01092021132420.fasta",
                  "./../results/EvAnSP_01092021133559/6_EvAnSP_01092021133559.fasta",
                  "./../results/EvAnSP_01092021134916/7_EvAnSP_01092021134916.fasta",
                  "./../results/EvAnSP_01092021140255/8_EvAnSP_01092021140255.fasta"]
    file_names2 = ["./../results/results/EvAnSP_01092021111859/2_EvAnSP_01092021111859.fasta",
                  "./../results/results/EvAnSP_01092021225442/3_EvAnSP_01092021225442.fasta",
                  "./../results/results/EvAnSP_01092021234619/4_EvAnSP_01092021234619.fasta",
                  "./../results/results/EvAnSP_01092021231100/5_EvAnSP_01092021231100.fasta",
                  "./../results/results/EvAnSP_02092021001206/6_EvAnSP_02092021001206.fasta",
                  "./../results/results/EvAnSP_01092021113444/7_EvAnSP_01092021113444.fasta",
                  "./../results/results/EvAnSP_01092021125805/8_EvAnSP_01092021125805.fasta"]
    o = []
    s = []
    for frame, file_name in enumerate(file_names2):
        output_list, same_output_list = get_evansp_output(file_name)
        # print(frame + 2)
        # print(min(output_list))
        # print(max(output_list))
        # print(np.mean(output_list))
        # print(np.std(output_list))
        # print(".")
        # print(min(same_output_list))
        # print(max(same_output_list))
        # print(np.mean(same_output_list))
        # print(np.std(same_output_list))
        # print("----")
        o.append(output_list)
        s.append((same_output_list))

    fig, (ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(7)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Probability')
    plt.ylabel('Number of LCRs')
    ax2.hist(x=s[0], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax2.set_ylabel("2")
    ax2.set_xticklabels([])
    ax3.hist(x=s[1], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax3.set_ylabel("3")
    ax3.set_xticklabels([])
    ax4.hist(x=s[2], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax4.set_ylabel("4")
    ax4.set_xticklabels([])
    ax5.hist(x=s[3], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax5.set_ylabel("5")
    ax5.set_xticklabels([])
    ax6.hist(x=s[4], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax6.set_ylabel("6")
    ax6.set_xticklabels([])
    ax7.hist(x=s[5], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax7.set_ylabel("7")
    ax7.set_xticklabels([])
    ax8.hist(x=s[6], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax8.set_ylabel("8")
    ax8.set_xticklabels([])

    # maxfreq = n.max()
    # plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig("./hist_2645_same_nolim.png")
    plt.show()

    o = []
    s = []
    for frame, file_name in enumerate(file_names):
        output_list, same_output_list = get_evansp_output(file_name)
        # print(frame + 2)
        # print(min(output_list))
        # print(max(output_list))
        # print(np.mean(output_list))
        # print(np.std(output_list))
        # print(".")
        # print(min(same_output_list))
        # print(max(same_output_list))
        # print(np.mean(same_output_list))
        # print(np.std(same_output_list))
        # print("----")
        o.append(output_list)
        s.append((same_output_list))

    fig, (ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(7)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Probability')
    plt.ylabel('Number of LCRs')
    ax2.hist(x=s[0], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax2.set_ylabel("2")
    ax2.set_xticklabels([])
    ax3.hist(x=s[1], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax3.set_ylabel("3")
    ax3.set_xticklabels([])
    ax4.hist(x=s[2], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax4.set_ylabel("4")
    ax4.set_xticklabels([])
    ax5.hist(x=s[3], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax5.set_ylabel("5")
    ax5.set_xticklabels([])
    ax6.hist(x=s[4], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax6.set_ylabel("6")
    ax6.set_xticklabels([])
    ax7.hist(x=s[5], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax7.set_ylabel("7")
    ax7.set_xticklabels([])
    ax8.hist(x=s[6], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax8.set_ylabel("8")
    ax8.set_xticklabels([])

    # maxfreq = n.max()
    # plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig("./hist_1000_same_nolim.png")
    plt.show()

    o = []
    s = []
    for frame, file_name in enumerate(file_names2):
        output_list, same_output_list = get_evansp_output(file_name)
        # print(frame + 2)
        # print(min(output_list))
        # print(max(output_list))
        # print(np.mean(output_list))
        # print(np.std(output_list))
        # print(".")
        # print(min(same_output_list))
        # print(max(same_output_list))
        # print(np.mean(same_output_list))
        # print(np.std(same_output_list))
        # print("----")
        o.append(output_list)
        s.append((same_output_list))

    fig, (ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(7)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Probability')
    plt.ylabel('Number of LCRs')
    ax2.hist(x=o[0], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax2.set_ylabel("2")
    ax2.set_xticklabels([])
    ax3.hist(x=o[1], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax3.set_ylabel("3")
    ax3.set_xticklabels([])
    ax4.hist(x=o[2], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax4.set_ylabel("4")
    ax4.set_xticklabels([])
    ax5.hist(x=o[3], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax5.set_ylabel("5")
    ax5.set_xticklabels([])
    ax6.hist(x=o[4], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax6.set_ylabel("6")
    ax6.set_xticklabels([])
    ax7.hist(x=o[5], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax7.set_ylabel("7")
    ax7.set_xticklabels([])
    ax8.hist(x=o[6], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax8.set_ylabel("8")
    ax8.set_xticklabels([])

    # maxfreq = n.max()
    # plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig("./hist_2645.png")
    plt.show()

    o = []
    s = []
    for frame, file_name in enumerate(file_names):
        output_list, same_output_list = get_evansp_output(file_name)
        # print(frame + 2)
        # print(min(output_list))
        # print(max(output_list))
        # print(np.mean(output_list))
        # print(np.std(output_list))
        # print(".")
        # print(min(same_output_list))
        # print(max(same_output_list))
        # print(np.mean(same_output_list))
        # print(np.std(same_output_list))
        # print("----")
        o.append(output_list)
        s.append((same_output_list))

    fig, (ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(7)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Probability')
    plt.ylabel('Number of LCRs')
    ax2.hist(x=o[0], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax2.set_ylabel("2")
    ax2.set_xticklabels([])
    ax3.hist(x=o[1], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax3.set_ylabel("3")
    ax3.set_xticklabels([])
    ax4.hist(x=o[2], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax4.set_ylabel("4")
    ax4.set_xticklabels([])
    ax5.hist(x=o[3], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax5.set_ylabel("5")
    ax5.set_xticklabels([])
    ax6.hist(x=o[4], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax6.set_ylabel("6")
    ax6.set_xticklabels([])
    ax7.hist(x=o[5], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax7.set_ylabel("7")
    ax7.set_xticklabels([])
    ax8.hist(x=o[6], bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    ax8.set_ylabel("8")
    ax8.set_xticklabels([])
    # maxfreq = n.max()
    # plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig("./hist_1000.png")
    plt.show()