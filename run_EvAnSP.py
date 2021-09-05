import click
import time

from src.parsers import Parser
from src.EvAnSP import EvAnSP


@click.command()
@click.option('--seq1', default='./LCR_sequences/examplary_LCRs.fasta',
              help='Path to file with first sequence.')
@click.option('--seq2', default='./LCR_sequences/examplary_LCRs.fasta',
              help='Path to file with second sequence.')
@click.option('--f', '--frame', default=5,
              help='Size of the frame used to analyse.')
@click.option('--all_frames', is_flag=True,
              help='Use this option, if you want to compere sequences' 
                   ' for all possible frames.')
@click.option('--mat', '--matrix',
              default='./prob_matrices/exp_prob_matrix',
              help='Matrix used to score alignments between sequences.')
def run(seq1, seq2, f, all_frames, mat):
    start_time = time.time()
    parser = Parser(seq1, seq2, mat)
    first_seq_dict, second_seq_dict, matrix = parser.get_parsered_input()
    print("Parser: --- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()

    program = EvAnSP(first_seq_dict, second_seq_dict, f, all_frames,  matrix)
    program.get_result()
    print("Evansp: --- %s seconds ---" % (time.time() - start_time))
    print("Task completed!")


if __name__ == "__main__":
    run()
