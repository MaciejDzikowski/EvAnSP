"""
Script by Joanna Ziemska-Legięcka.
Fixed and refactored by Maciej Dzikowski.
"""

import glob
import os
import align_mafft as am
import get_lcrs_stats as stats

from typing import Callable, Dict


def load_sequences(path: str, header_type: Callable) -> Dict[str, Dict[str, str]]:
    # returns: {family_name:group_num: {header: seq}
    groups = {}
    for filename in glob.glob(os.path.join(path, '*.fasta')):
        with open(os.path.join(os.getcwd(), filename), 'r') as f:
            seqs = {}
            header = False
            counter = 0
            for line in f:
                if line.startswith(">"):
                    header = header_type(line)
                    seqs[header.__str__()] = ""
                elif line.startswith("-----"):
                    groups[header.get_family_name() + ":" + str(counter)] = seqs
                    seqs = {}
                    header = False
                    counter += 1
                elif line:
                    if header:
                        seqs[header.__str__()] += line.strip()
    return groups


def divide_cluster_by_aminoacid_composition_2(groups: Dict[str, Dict[str, str]], value: float) -> Dict[str, Dict[str, str]]:
    for cluster_name, cluster_seq in list(groups.items()):
        divide_one_cluster_by_aminoacid_composition_2(value, cluster_name, cluster_seq, groups)
    return groups


def count_aminoacid_composition(sequence: str) -> Dict[str, float]:
    return {aa: sequence.count(aa) / len(sequence) for aa in set(sequence)}


def count_diffrence(hist_sequence1: Dict[str, float], hist_sequence2: Dict[str, float]) -> float:
    diffrence_counter = 0
    # hist_sequence1 {"A":0.1, "G":0.9}

    # hist_sequence2 {S:0,1, A:0.1, H:0,8}
    for i in list(list(hist_sequence1.keys()) + list(hist_sequence2.keys())):
        if i not in hist_sequence1:
            diffrence_counter += hist_sequence2[i]
        elif i not in hist_sequence2:
            diffrence_counter += hist_sequence1[i]
        elif i in hist_sequence2 and i in hist_sequence1:
            if hist_sequence1[i] > hist_sequence2[i]:
                diffrence_counter += hist_sequence1[i] - hist_sequence2[i]
            else:
                diffrence_counter += hist_sequence2[i] - hist_sequence1[i]
    return 100 * diffrence_counter / 2


def divide_one_cluster_by_aminoacid_composition_2(value: float, cluster_name: str, cluster_seqs: Dict[str, str],
                                                  groups: Dict[str, Dict[str, str]]):
    new_cluster_proteins = {}
    if len(cluster_seqs) > 1:
        centroid = cluster_seqs[next(iter(cluster_seqs))]

        centroid_composition = count_aminoacid_composition(centroid)
        for seq_header, seq_sequence in cluster_seqs.items():
            sequence_composition = count_aminoacid_composition(seq_sequence)
            diffrence_value = count_diffrence(sequence_composition, centroid_composition)
            if diffrence_value is not None and diffrence_value > value:
                new_cluster_proteins[seq_header] = seq_sequence
                groups[cluster_name] = {i: j for i, j in groups[cluster_name].items() if i != seq_header}

        if new_cluster_proteins:
            if "|" in cluster_name:
                cluster_name, cluster_number = cluster_name.rsplit("|", 1)
                cluster_number = int(cluster_number) + 1
            else:
                cluster_name, cluster_number = cluster_name, 0
            cluster_value = f"{cluster_name}|{cluster_number}"

            if cluster_value not in groups:
                groups[cluster_value] = new_cluster_proteins
            else:
                cluster_value = f"{cluster_name}|{cluster_number}_|0"
                groups[cluster_value] = new_cluster_proteins
            divide_one_cluster_by_aminoacid_composition_2(value, cluster_value, new_cluster_proteins, groups)


def save_maffted_clusters(clusters: Dict[str, Dict[str, str]], output_name: str):
    for cluster in clusters:
        temp_file_name = am.write_temp_file(clusters[cluster], output_name)
        mafft_file_name = am.run_mafft(temp_file_name, cluster, output_name)
        os.remove(temp_file_name)

        mafft_output = am.get_mafft_output(mafft_file_name, stats.PfamHeader)
        os.remove(mafft_file_name)

        test_len = len(list(mafft_output.values())[0])

        cluster_file_name = "./{}/{}.txt".format(output_name, cluster)
        if os.path.dirname(cluster_file_name):
            os.makedirs(os.path.dirname(temp_file_name), exist_ok=True)
        file = open(cluster_file_name, "w")
        for header, seq in mafft_output.items():
            if len(seq) != test_len:
                raise Exception('Maffted sequences has no equal length!')
            file.write(">" + header.__str__() + "\t")
            file.write(seq.upper() + "\n")
        file.close()


if __name__ == '__main__':
    value = 60
    output_name = "Pfam_LCRs_clustered"

    groups = load_sequences('./test/', stats.PfamHeaderLCR)  # mafft_Pfam_LCRs_grouped
    clusters = divide_cluster_by_aminoacid_composition_2(groups, value)
    save_maffted_clusters(clusters, output_name)
