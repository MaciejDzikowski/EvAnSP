import os
import scripts.get_lcrs_stats as stats

from datetime import datetime
from src.parsers import Sequence
from typing import Dict, List


class EvAnSP:
    def __init__(self, first_seq_dict: List[Sequence], second_seq_dict: List[Sequence],
                 frame_size: int, all_frames: bool, matrix: Dict[str, float]):
        self.first_seq_dict = first_seq_dict
        self.second_seq_dict = second_seq_dict
        self.frame_size = frame_size
        self.all_frames = all_frames
        self.matrix = matrix

    def get_probability(self, seq1: str, seq2: str, frame_size: int) -> float:
        score = 0
        matrix = self.matrix
        for i in range(len(seq1) - frame_size + 1):
            for j in range(len(seq2) - frame_size + 1):
                frame_score = 0
                for aa in range(frame_size):
                    a_in_a_prob = matrix[seq1[i + aa] + seq2[j + aa]]
                    frame_score += a_in_a_prob
                score += frame_score / frame_size
        probability = score / ((len(seq1) - frame_size + 1) * (len(seq2) - frame_size + 1))
        return probability

    def get_result(self):
        first_seq_dict = self.first_seq_dict
        second_seq_dict = self.second_seq_dict
        all_frames = self.all_frames
        dir_name = datetime.now().strftime("EvAnSP_%d%m%Y%H%M%S")

        frames = [i for i in range(1, min(int(len(max(first_seq_dict, key=lambda k: len(k.seq)).seq)),
                                          int(len(max(second_seq_dict, key=lambda k: len(k.seq)).seq))))] \
            if all_frames else [self.frame_size]
        print("We will use {} frames.".format(len(frames)))
        # print("Comparing sequences for frame:", frame_size)
        results_list = {}
        num_first_seqs = len(first_seq_dict)
        for i, seq1 in enumerate(first_seq_dict):
            # print("[{}/{}] {}".format(i + 1, num_first_seqs, seq1))
            for j, seq2 in enumerate(second_seq_dict):
                for frame_size in frames:
                    if frame_size <= min(len(seq1.seq), len(seq2.seq)):
                        probability = self.get_probability(seq1.seq, seq2.seq, frame_size)
                        if frame_size not in results_list:
                            results_list[frame_size] = {}
                        if seq1 not in results_list[frame_size]:
                            results_list[frame_size][seq1] = {}
                        if seq2 not in results_list[frame_size][seq1]:
                            results_list[frame_size][seq1][seq2] = {}
                        results_list[frame_size][seq1][seq2] = probability
        for frame in results_list:
            save_results(results_list[frame], frame, dir_name)


def save_results(results_list: dict, frame_size: int, dir_name: str):
    print("Saving results for frame: {}".format(frame_size))
    file_name = datetime.now().strftime("./results/{}/{}_{}.fasta".format(dir_name, str(frame_size), dir_name))
    if os.path.dirname(file_name):
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
    file = open(file_name, "w")
    for seq1 in results_list:
        for seq2 in results_list[seq1]:
            file.write(seq1.header.__str__() + "\n")
            file.write(seq2.header.__str__() + "\n")
            file.write(seq1.seq + " ~ " + seq2.seq + "\n")
            file.write("probability: " + str(results_list[seq1][seq2]) + " frame size: "
                       + str(frame_size) + "\n")
            file.write("\n")
    file.close()
