"""
Script by Maciej Dzikowski.
"""

import click
import os
import get_lcrs_sequences as gls

from datetime import datetime


def get_swiss_prot_seqs(file: str) -> dict:
    return {line.strip().strip(">").split()[0].split("|")[-1]: line.strip().strip(">")
            for line in open(file) if line.startswith(">")}


def get_common_pfam_seqs(file: str, seqs: dict) -> dict:
    common_seq = {}
    checker = False
    with open(file) as f:
        for line in f:
            if line.startswith(">"):
                if line.strip().strip(">").split("/")[0] in seqs:
                    checker = (line.strip(';\n') + " OS="
                               + seqs[line.strip().strip(">").split("/")[0]].split("OS=")[1].split(" LCR:begin=")[0].split(" GN=")[0])
                    print(checker)
                else:
                    checker = False
            elif checker:
                if checker in common_seq.keys():
                    common_seq[checker] += line.strip()
                else:
                    common_seq[checker] = line.strip()
    return common_seq


def save_new_file_from_dict(seq_dict: dict, file_name: str):
    if os.path.dirname(file_name):
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
    with open(file_name, "w") as file:
        for record in seq_dict:
            file.write(record + "\n")
            file.write(seq_dict[record] + "\n")


@click.command()
@click.option("--sp", type=str, default="./LCR_sequences/swiss_prot_LCRs.fasta",
              help="Path to Swiss-Prot .fasta database with LCRs.")
@click.option("--pfam", type=str, default="./sequences/Pfam-A.fasta", help="Path to Pfam .fasta database.")
@click.option("--o", type=str, default=datetime.now().strftime("./LCR_sequences/Pfam-SP_LCRs_%d%m%Y%H%M%S.fasta"),
              help="Output file name.")
@click.option("--seg_path", type=str, default="segmasker", help="Path to segmasker.")
@click.option("--seg_param", type=str, default="{'locut': 1.5, 'hicut': 1.8, 'window': 15}",
              help="Parameters of seg in dictionary. For example: \"{'locut': 1.5, 'hicut': 1.8, 'window': 15}\"")
@click.pass_context
def run(ctx, sp: str, pfam: str, o: str, seg_path: str, seg_param: dict):
    print('[1/2] Finding common sequences...')
    try:
        temp_file = "./temporary_common_seqs_file.fasta"
        save_new_file_from_dict(get_common_pfam_seqs(pfam, get_swiss_prot_seqs(sp)), temp_file)
    except Exception as e:
        print('Finding common sequences: failed. Error:\n', e)
    else:
        print('[2/2] Running SEG...')
        try:
            ctx.invoke(gls.run, sequence_db=temp_file, output_file=o, seg_path=seg_path, seg_param=seg_param)
        except Exception as e:
            print('Running SEG: failed. Error:\n', e)
        else:
            print("Task completed!")
        finally:
            os.remove(temp_file)


if __name__ == "__main__":
    run()
