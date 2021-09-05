codons = {
    "UUU": "F",
    "CUU": "L",
    "AUU": "I",
    "GUU": "V",
    "UUC": "F",
    "CUC": "L",
    "AUC": "I",
    "GUC": "V",
    "UUA": "L",
    "CUA": "L",
    "AUA": "I",
    "GUA": "V",
    "UUG": "L",
    "CUG": "L",
    "AUG": "M",
    "GUG": "V",
    "UCU": "S",
    "CCU": "P",
    "ACU": "T",
    "GCU": "A",
    "UCC": "S",
    "CCC": "P",
    "ACC": "T",
    "GCC": "A",
    "UCA": "S",
    "CCA": "P",
    "ACA": "T",
    "GCA": "A",
    "UCG": "S",
    "CCG": "P",
    "ACG": "T",
    "GCG": "A",
    "UAU": "Y",
    "CAU": "H",
    "AAU": "N",
    "GAU": "D",
    "UAC": "Y",
    "CAC": "H",
    "AAC": "N",
    "GAC": "D",
    "UAA": "Stop",
    "CAA": "Q",
    "AAA": "K",
    "GAA": "E",
    "UAG": "Stop",
    "CAG": "Q",
    "AAG": "K",
    "GAG": "E",
    "UGU": "C",
    "CGU": "R",
    "AGU": "S",
    "GGU": "G",
    "UGC": "C",
    "CGC": "R",
    "AGC": "S",
    "GGC": "G",
    "UGA": "Stop",
    "CGA": "R",
    "AGA": "R",
    "GGA": "G",
    "UGG": "W",
    "CGG": "R",
    "AGG": "R",
    "GGG": "G"}
number_of_codons = 64
nucleotides = {"U": 0.22, "A": 0.303, "C": 0.217, "G": 0.261}  # Why isn't sum equal 1?


def get_aas_prob(codons, nucleotides, number_of_codons):
    aas_probs = {}
    for codon in codons:
        sum = 0
        for nucleotide in codon:
            sum += nucleotides[nucleotide]
        prob = sum / number_of_codons
        if codons[codon] in aas_probs:
            aas_probs[codons[codon]] += prob
        else:
            aas_probs[codons[codon]] = prob

    X_prob = 0
    for aa in aas_probs:
        X_prob += aas_probs[aa]
    aas_probs["X"] = X_prob / len(aas_probs)
    aas_probs["B"] = aas_probs["N"] + aas_probs["D"]
    aas_probs["Z"] = aas_probs["E"] + aas_probs["Q"]
    return aas_probs


print(get_aas_prob(codons, nucleotides, number_of_codons))
