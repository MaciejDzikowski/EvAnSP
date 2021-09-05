import click

from math import log2, log10
from datetime import datetime


class Creator:
    def __init__(self, db_file, aas_replacing_prob, file_name):
        self.db_file = db_file
        self.aas_replacing_prob = aas_replacing_prob
        self.file_name = file_name

        self.aas_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L",
                         "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z",
                         "x"]

        self.parsed_data = self.parser()
        self.counted_aas, self.num_aas = self.count_aas()
        self.matrix = self.create_matrix()

    def parser(self):
        parsed_data = []
        for line in open(self.db_file):
            if not line.startswith(">"):
                parsed_data.append(line.strip().split())
        self.db_file.close()
        return parsed_data

    def count_aas(self):
        counted_aas = {}
        num_aas = 0
        for sequence in self.parsed_data:
            for nucleotide in sequence:
                if nucleotide in counted_aas:
                    counted_aas[nucleotide] += 1
                else:
                    counted_aas[nucleotide] = 1
                num_aas += 1
        return counted_aas, num_aas

    def create_matrix(self):
        matrix = []
        la = (1 / 2) * log10(2)
        for i, aa1 in enumerate(self.aas_list):
            scores = []
            for j, aa2 in enumerate(self.aas_list):
                s = (1 / la) * log10(2) * log2(self.aas_replacing_prob[i][j]
                                               / ((aa1 / self.num_aas)
                                                  * (aa2 / self.num_aas)))
                scores.append(s)
            scores.append(min(scores))  # adds value for "*"
            matrix.append(scores)
        return matrix

    def save_matrix(self):
        matrix = self.matrix
        # try:
        file_name = datetime.now().strftime(
            "./prob_matrices/" + self.file_name)
        file = open(file_name, "w")
        file.write("# Matrix made by matrix_creator from EvAnSP" + "\n")
        file.write("# on " + datetime.now().strftime("%d/%m/%Y %H:%M.") + "\n")
        file.write("# BLOSUM Clustered Scoring Matrix" + "\n")
        file.write("# Used database: " + self.db_file + "\n")
        file.write("# * column uses minimum score" + "\n")

        file.write(" ")
        for aa_abb in self.aas_list:
            file.write("  " + aa_abb)
        file.write("\n")

        for i, aa in enumerate(matrix):
            file.write(self.aas_list[i])
            for score in matrix[i]:
                file.write(format(score, ' 3d'))
            file.write("\n")
        file.close()


@click.command()
@click.option('--db_file', default='./LCR_sequences/swiss_prot_LCRs.fasta',
              help='Sequences you want to use to create a scoring matrix.')
@click.option('--arp', default='',
              help='Probability of two amino acids replacing each other.'
                   'Order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P,'
                   ' S, T, W, Y, V, B, Z, X')
@click.option('--file_name', '--name',
              default=datetime.now().strftime("/smatrix_%d%m%Y%H%M%S"),
              help='Use this option, if you want to name your output file.')
def run(db_file, arp, file_name):
    creator = Creator(db_file, arp, file_name)
    creator.save_matrix()


if __name__ == "__main__":
    run()
