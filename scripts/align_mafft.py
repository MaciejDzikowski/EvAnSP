import get_lcrs_stats as stats
import numpy as np
import os

from typing import Callable, Dict, Tuple


def get_full_seqs_per_family(seqs: Dict[stats.PfamHeader, str]) -> Dict[str, Dict[stats.PfamHeader, str]]:
    families = {}
    for seq in seqs:
        family = seq.get_family_name()
        if family in families:
            families[family][seq] = seqs[seq]
        else:
            families[family] = {seq: seqs[seq]}
    return families


def get_common_pfam_seqs_per_family(file: str, families_seqs: Dict[str, Dict[stats.PfamHeader, str]]) -> Dict[stats.PfamHeader, str]:
    common_seq = {}
    seqs_names = set([seq.get_identification_str() for family in families_seqs for seq in families_seqs[family]])
    checker = False
    with open(file) as f:
        for line in f:
            if line.startswith(">"):
                if line.strip().strip(">").split()[0] in seqs_names:
                    checker = stats.PfamHeader(line.strip(';\n'))
                    print(checker)
                else:
                    checker = False
            elif checker:
                if checker in common_seq.keys():
                    common_seq[checker] += line.strip()
                else:
                    common_seq[checker] = line.strip()
    return common_seq


def write_temp_file(familt_prot_seq: dict, output_name: str) -> str:
    temp_file_name = "./temp_for_{}.fasta".format(output_name)
    if os.path.dirname(temp_file_name):
        os.makedirs(os.path.dirname(temp_file_name), exist_ok=True)
    file = open(temp_file_name, "w")
    for header, seq in familt_prot_seq.items():
        file.write(">" + header.__str__() + "\n")
        file.write(seq + "\n")
    return temp_file_name


def run_mafft(file_name: str, family_name: str, output_name: str) -> str:
    output_name = "./{}/mafft_{}. fasta".format(output_name, family_name)
    if os.path.dirname(output_name):
        os.makedirs(os.path.dirname(output_name), exist_ok=True)
    command = 'mafft --auto --inputorder "{}" > "{}"'.format(file_name, output_name)
    os.system(command)
    if os.stat(output_name).st_size == 0:
        os.remove(output_name)
        command = 'mafft --auto --inputorder --anysymbol "{}" > "{}"'.format(file_name, output_name)
        os.system(command)
    return output_name


def get_mafft_output(file_name: str, header_type: Callable) -> Dict[stats.Header, str]:
    return stats.get_seqs_dict(file_name, header_type)


def mark_lcrs(mafft_output: Dict[stats.Header, str], lcrs: Dict[stats.PfamHeaderLCR, str]) -> Dict[str, str]:
    marked_html = {}
    for header in mafft_output:
        seq = mafft_output[header]
        slicing_pos = []
        for lcr in lcrs:
            if lcr.get_entry_name() == header.get_entry_name():
                slicing_pos.append(lcr.get_lcr_start())
                slicing_pos.append(lcr.get_lcr_end() + 1)
        slicing_pos = np.array(slicing_pos)
        marked_seq = ""
        for i, aa in enumerate(seq):
            if i == 0 and 0 in slicing_pos and aa != "-":
                marked_seq += "$"
            if aa == "-":
                slicing_pos += 1
            marked_seq += aa
            if i + 1 in slicing_pos and i + 1 != len(seq) and seq[i+1] != "-":
                marked_seq += "$"
        marked_html[header.__str__()] = marked_seq
    return marked_html


def write_html_file(parsed_mafft_output: Dict[str, str], family_name: str,
                    output_name: str) -> Tuple[str, str]:
    if parsed_mafft_output:
        file_name = "./{}/mafft_{}.html".format(output_name, family_name)
        if os.path.dirname(file_name):
            os.makedirs(os.path.dirname(file_name), exist_ok=True)
        file = open(file_name, "w")
        file.write("<header>" + "\n")
        file.write("<h3><a style=\"color: black;\" href=\"./../{}.html\">Choose another family</a></h5>".format(output_name) + "\n")
        file.write("</header>" + "\n")
        file.write("<hr>" + "\n")
        file.write("<table>" + "\n")
        file.write("\t" + "<tr>" + "\n")
        file.write("\t" + "\t" + "<th style=\"text-align: left; background-color: gainsboro;\">" + "<nobr>" + "Header" + "<nobr>" + "</th>" + "\n")
        file.write("\t" + "\t" + "<th style=\"text-align: left; background-color: gainsboro;\">" + "Sequence" + "</th>" + "\n")
        file.write("\t" + "</tr>" + "\n")

        is_even = False
        for header, seq in parsed_mafft_output.items():
            if is_even:
                file.write("\t" + "<tr>" + "\n")
            else:
                file.write("\t" + "<tr style=\" background-color: whitesmoke;\">" + "\n")
            file.write("\t" + "\t" + "<td>" + "<nobr>" + header + "</td>" + "\n")
            file.write("\t" + "\t" + "<td>" + "<nobr>" + "\n")
            is_with_tag = False
            for aa in seq:
                if aa == "$":
                    is_with_tag = not is_with_tag
                    continue
                if is_with_tag and aa != "-":
                    file.write("\t" + "\t" + "\t" + "<span style=\"text-align: center; display: inline-block; width: 15px; color: red;\">" + aa + "</span>" + "\n")
                else:
                    file.write("\t" + "\t" + "\t" + "<span style=\"text-align: center; display: inline-block; width: 15px;\">" + aa + "</span>" + "\n")
            file.write("\t" + "\t" + "</td>" + "\n")
            file.write("\t" + "</nobr>" + "</tr>" + "\n")
            is_even = not is_even
        file.write("</table>" + "\n")
        file.close()
    else:
        print(family_name)
        file_name = "./mafft/aErr_{}.html".format(family_name)
        if os.path.dirname(file_name):
            os.makedirs(os.path.dirname(file_name), exist_ok=True)
        file = open(file_name, "w")
        file.close()
    return file_name, family_name


def create_index_html(htmls_families: list, output_name: str):
    file_name = "./{}.html".format(output_name)
    if os.path.dirname(file_name):
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
    file = open(file_name, "w")
    file.write("<header>" + "\n")
    file.write("<h3 style=\"color: black; font-weight: bold; text-align: center;\">Choose family to see an alignment:</h5>" + "\n")
    file.write("</header>" + "\n")
    file.write("<hr>" + "\n")
    file.write("<head>" + "\n")
    file.write("\t" + "<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css\" integrity=\"sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm\" crossorigin=\"anonymous\"/>" + "\n")
    file.write("</head>" + "\n")
    file.write("<body><divclass=\"container-fluid\">" + "\n")
    for i, t in enumerate(htmls_families):
        if i % 3 == 0:
            file.write("\t" + "<div class=\"row\">" + "\n")
        file.write("\t" + "\t" + "<div class=\"col-lg-4 col-md-6 col-xs-12\">" + "\n")
        file.write("\t" + "\t" + "\t" + "<li><a style=\"color: black; text-align: center;\" href=\"{}\">{}</a></li>".format(t[0], t[1]) + "\n")
        file.write("\t" + "\t" + "</div>" + "\n")
        if i % 3 == 2:
            file.write("\t" + "</div>" + "\n")
    file.write("</div></body>" + "\n")
    file.close()


if __name__ == "__main__":
    # lcrs_file_name = "./../LCR_sequences/Pfam_LCRs.fasta"
    lcrs_file_name = "./../LCR_sequences/test_mafft_lcrs.fasta"
    header_type = stats.PfamHeaderLCR
    # pfam_file_name = "./../sequences/Pfam-A.fasta"
    pfam_file_name = "./../LCR_sequences/test_mafft.fasta"
    output_name = "mafft_Pfam_LCRs_new"
    lcr_families = get_full_seqs_per_family(stats.get_seqs_dict(lcrs_file_name, header_type))
    mafft_lcr_families = {k: v for k, v in lcr_families.items()
                          if len(v) > 1 and len(set([i.get_identification_str() for i in v])) > 1}
    prot_seqs = get_common_pfam_seqs_per_family(pfam_file_name, mafft_lcr_families)

    htmls_families = []
    for family in mafft_lcr_families:
        family_prot_seqs = {header: prot_seqs[header] for header in prot_seqs
                            if header.__str__() in [i.__str__().split(" LCR")[0]
                                                    for i in mafft_lcr_families[family]]}
        temp_file_name = write_temp_file(family_prot_seqs, output_name)
        mafft_file_name = run_mafft(temp_file_name, family, output_name)
        os.remove(temp_file_name)

        htmls_families.append(
            write_html_file(
                mark_lcrs(
                    get_mafft_output(mafft_file_name, stats.PfamHeader),
                    mafft_lcr_families[family]
                ),
                family,
                output_name
            )
        )
        os.remove(mafft_file_name)
    create_index_html(htmls_families, output_name)

## TODO
# use click module
# add prev/ next family buttons