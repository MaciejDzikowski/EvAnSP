"""
Script by Joanna Ziemska-Legięcka.
Improved and refactored by Maciej Dzikowski.
"""

import click
import os

from subprocess import Popen, PIPE


class LCRsSequences:
    def __init__(self, input_file: str, output_file: str,
                 seg_parameters: dict=None, seg_path: str="segmasker"):
        self.input = input_file
        self.output_file = output_file
        self.seg_path = seg_path
        self.seg_parameters = seg_parameters

        self.output = None
        self.proteins = {}

    def prepare_params(self) -> str:
        txt = ""
        if self.seg_parameters is None:
            txt = "-locut 1.5 -hicut 1.8 -window 15"
        else:
            for value, key in self.seg_parameters.items():
                txt += " -" + value + " " + str(key)
        return txt

    def run_seg(self):
        input_fasta = open(self.input, "r").read()
        params = self.seg_path + " " + self.prepare_params()
        p = Popen(params.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate(input=input_fasta.encode("ascii"))
        if stderr:
            raise Exception('SegmaskerError:\n' + stderr.decode())
        self.output = stdout.decode()
        parsed_output = self.parse_output()
        self.write_output(parsed_output)

    def parse_output(self) -> dict:
        retval = {}
        protein_list = self.read_input()
        for line in self.output.splitlines():
            if line.startswith(">"):
                seq = protein_list[line]
                retval[line] = []
                last_header = line
            else:
                line_items = line.split("-")
                beg = int(line_items[0].strip())
                end = int(line_items[1].strip())
                region = seq[beg:end+1]
                if last_header in retval.keys():
                    retval[last_header].append({"region": region, "end": end, "beg": beg})
                else:
                    retval[last_header] = [{"region": region, "end": end, "beg": beg}]
        return retval

    def read_input(self) -> dict:
        res = {}
        header = False
        with open(self.input) as f:
            for line in f:
                if line.startswith(">"):
                    header = line.strip().strip(";")
                    res[header] = ""
                elif line:
                    if header:
                        res[header] += line.strip()
        return res

    def write_output(self, parsed_output: dict):
        if os.path.dirname(self.output_file):
            os.makedirs(os.path.dirname(self.output_file), exist_ok=True)
        file = open(self.output_file, "w")
        for header, regions in parsed_output.items():
            if len(regions) > 0:
                for region in regions:
                    file.write(header + " LCR:begin=" + str(region['beg']) + ", end=" + str(region['end']) + "\n")
                    file.write(region['region'] + "\n")
        file.close()


@click.command()
@click.option("--sequence_db", type=str, default="./examples/LCRs.fasta", help="Path to fasta database with LCRs.")
@click.option("--output_file", type=str, default="./seg_tmp.fasta", help="Output with clusters of LCRs.")
@click.option("--seg_path", type=str, default="segmaskr", help="Path to segmasker.")
@click.option("--seg_param", type=str, default="{'locut': 1.5, 'hicut': 1.8, 'window': 15}",
              help="Parameters of seg in dictionary. For example: \"{'locut': 1.5, 'hicut': 1.8, 'window': 15}\"")
def run(sequence_db: str, output_file: str, seg_param: str, seg_path: str):
    if seg_param is not None:
        tmp = LCRsSequences(sequence_db, output_file, eval(seg_param), seg_path)
    tmp.run_seg()


if __name__ == "__main__":
    run()
