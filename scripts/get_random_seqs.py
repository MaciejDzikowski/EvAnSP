import os
import random
import src.parsers as p


if __name__ == "__main__":
    # draw = 1000
    # file_name = "./../LCR_sequences/swiss_prot_LCRs_1000.fasta"

    lcrs = p.create_seq_list("./../LCR_sequences/swiss_prot_LCRs_052021.fasta")
    lcrs_num = len(lcrs)

    draw = int(lcrs_num / 10)
    file_name = "./../LCR_sequences/swiss_prot_LCRs_{}.fasta".format(draw)

    print(draw)

    random_indexes = []
    while len(random_indexes) != draw:
        random_int = random.randint(0, lcrs_num-1)
        if random_int not in random_indexes:
            random_indexes.append(random_int)

    if os.path.dirname(file_name):
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
    with open(file_name, "w") as file:
        for i in random_indexes:
            file.write(">" + lcrs[i].__str__() + "\n")
            file.write(lcrs[i].seq + "\n")
