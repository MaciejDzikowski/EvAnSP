import scripts.get_lcrs_stats as st


# Test header classes
eg_uniprot_headers = [
    ">sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1",
    ">sp|Q6GZW6|009L_FRG3G Putative helicase 009L OS=Frog virus 3 (isolate Goorha) OX=654924 PE=4 SV=1"
]

result_uniprot_headers = [
    ("sp", "Q6GZX4", "001R_FRG3G", "Putative transcription factor 001R", "Frog virus 3 (isolate Goorha)", "654924", "FV3-001R", 4, 1),
    ("sp", "Q6GZW6", "009L_FRG3G", "Putative helicase 009L", "Frog virus 3 (isolate Goorha)", "654924", None, 4, 1)
]

eg_pfam_headers = [
    ">A0A089UQM1_9ENTR/154-186 A0A089UQM1.1 PF10417.7;1-cysPrx_C;",
    ">A0A0D2WKE7_CAPO3/1080-1115 A0A0D2WKE7.1 PF10417.7;1-cysPrx_C;",
    ">G3RSJ9_GORGO/210-245 G3RSJ9.1 PF10417.7;1-cysPrx_C;"
]

result_pfam_headers = [
    ("A0A089UQM1_9ENTR", 154, 186, "A0A089UQM1.1", "PF10417", "1-cysPrx_C"),
    ("A0A0D2WKE7_CAPO3", 1080, 1115, "A0A0D2WKE7.1", "PF10417", "1-cysPrx_C"),
    ("G3RSJ9_GORGO", 210, 245, "G3RSJ9.1", "PF10417", "1-cysPrx_C")
]

eg_uniprot_lcr_headers = [
    ">sp|O55724|106L_IIV6 Uncharacterized protein 106L OS=Invertebrate iridescent virus 6 OX=176652 GN=IIV6-106L PE=4 SV=1 LCR:begin=24, end=56",
    ">sp|P39413|AEF1_DROME Adult enhancer factor 1 OS=Drosophila melanogaster OX=7227 PE=2 SV=1 LCR:begin=107, end=130"
]

result_uniprot_lcr_headers = [
    ("sp", "O55724", "106L_IIV6", "Uncharacterized protein 106L", "Invertebrate iridescent virus 6", "176652", "IIV6-106L", 4, 1, 24, 56),
    ("sp", "P39413", "AEF1_DROME", "Adult enhancer factor 1", "Drosophila melanogaster", "7227", None, 2, 1, 107, 130)
]

eg_pfam_lcr_headers = [
    ">A0A0G2FN91_9PEZI/248-391 A0A0G2FN91.1 PF03171.18;2OG-FeII_Oxy LCR:begin=43, end=58",
    ">A0A0F8WC61_9EURO/142-401 A0A0F8WC61.1 PF13532.4;2OG-FeII_Oxy_2 LCR:begin=176, end=191"
]

result_pfam_lcr_headers = [
    ("A0A0G2FN91_9PEZI", 248, 391, "A0A0G2FN91.1", "PF03171", "2OG-FeII_Oxy", 43, 58),
    ("A0A0F8WC61_9EURO", 142, 401, "A0A0F8WC61.1", "PF13532", "2OG-FeII_Oxy_2", 176, 191)
]

eg_pfam_in_sp_lcr_headers = [
    ">RSSA_CULQU/200-285 B0X6V0.2 PF16122.3;40S_SA_C OS=Culex quinquefasciatus OX=7176 LCR:begin=41, end=62",
    ">BFR1_SCHPO/20-158 P41820.1 PF14510.4;ABC_trans_N OS=Schizosaccharomyces pombe (strain 972 / ATCC 24843) OX=284812 LCR:begin=62, end=75"
]

result_pfam_in_sp_lcr_headers = [
    ("RSSA_CULQU", 200, 285, "B0X6V0.2", "PF16122", "40S_SA_C", "Culex quinquefasciatus", "7176", 41, 62),
    ("BFR1_SCHPO", 20, 158, "P41820.1", "PF14510", "ABC_trans_N", "Schizosaccharomyces pombe (strain 972 / ATCC 24843)", "284812", 62, 75)
]

print('Uniprot headers...')
for i, header in enumerate(eg_uniprot_headers):
    tested_header = st.UniprotHeader(header)
    if tested_header.get_database_id() != result_uniprot_headers[i][0]:
        print(tested_header.get_database_id(), "!=", result_uniprot_headers[i][0])
    if tested_header.get_unique_identifer() != result_uniprot_headers[i][1]:
        print(tested_header.get_unique_identifer(), "!=", result_uniprot_headers[i][1])
    if tested_header.get_entry_name() != result_uniprot_headers[i][2]:
        print(tested_header.get_entry_name(), "!=", result_uniprot_headers[i][2])
    if tested_header.get_protein_name() != result_uniprot_headers[i][3]:
        print(tested_header.get_protein_name(), "!=", result_uniprot_headers[i][3])
    if tested_header.get_organism_name() != result_uniprot_headers[i][4]:
        print(tested_header.get_organism_name(), "!=", result_uniprot_headers[i][4])
    if tested_header.get_organism_identifier() != result_uniprot_headers[i][5]:
        print(tested_header.get_organism_identifier(), "!=", result_uniprot_headers[i][5])
    if tested_header.get_gene_name() != result_uniprot_headers[i][6]:
        print(tested_header.get_gene_name(), "!=", result_uniprot_headers[i][6])
    if tested_header.get_protein_existence() != result_uniprot_headers[i][7]:
        print(tested_header.get_protein_existence(), "!=", result_uniprot_headers[i][7])
    if tested_header.get_sequence_version() != result_uniprot_headers[i][8]:
        print(tested_header.get_sequence_version(), "!=", result_uniprot_headers[i][8])

print('Pfam headers...')
for i, header in enumerate(eg_pfam_headers):
    tested_header = st.PfamHeader(header)
    if tested_header.get_entry_name() != result_pfam_headers[i][0]:
        print(tested_header.get_entry_name(), "!=", result_pfam_headers[i][0])
    if tested_header.get_prot_start() != result_pfam_headers[i][1]:
        print(tested_header.get_prot_start(), "!=", result_pfam_headers[i][1])
    if tested_header.get_prot_end() != result_pfam_headers[i][2]:
        print(tested_header.get_prot_end(), "!=", result_pfam_headers[i][2])
    if tested_header.get_unique_identifer() != result_pfam_headers[i][3]:
        print(tested_header.get_unique_identifer(), "!=", result_pfam_headers[i][3])
    if tested_header.get_family_id() != result_pfam_headers[i][4]:
        print(tested_header.get_family_id(), "!=", result_pfam_headers[i][4])
    if tested_header.get_family_name() != result_pfam_headers[i][5]:
        print(tested_header.get_family_name(), "!=", result_pfam_headers[i][5])

print('Uniprot LCR headers...')
for i, header in enumerate(eg_uniprot_lcr_headers):
    tested_header = st.UniprotHeaderLCR(header)
    if tested_header.get_database_id() != result_uniprot_lcr_headers[i][0]:
        print(tested_header.get_database_id(), "!=", result_uniprot_lcr_headers[i][0])
    if tested_header.get_unique_identifer() != result_uniprot_lcr_headers[i][1]:
        print(tested_header.get_unique_identifer(), "!=", result_uniprot_lcr_headers[i][1])
    if tested_header.get_entry_name() != result_uniprot_lcr_headers[i][2]:
        print(tested_header.get_entry_name(), "!=", result_uniprot_lcr_headers[i][2])
    if tested_header.get_protein_name() != result_uniprot_lcr_headers[i][3]:
        print(tested_header.get_protein_name(), "!=", result_uniprot_lcr_headers[i][3])
    if tested_header.get_organism_name() != result_uniprot_lcr_headers[i][4]:
        print(tested_header.get_organism_name(), "!=", result_uniprot_lcr_headers[i][4])
    if tested_header.get_organism_identifier() != result_uniprot_lcr_headers[i][5]:
        print(tested_header.get_organism_identifier(), "!=", result_uniprot_lcr_headers[i][5])
    if tested_header.get_gene_name() != result_uniprot_lcr_headers[i][6]:
        print(tested_header.get_gene_name(), "!=", result_uniprot_lcr_headers[i][6])
    if tested_header.get_protein_existence() != result_uniprot_lcr_headers[i][7]:
        print(tested_header.get_protein_existence(), "!=", result_uniprot_lcr_headers[i][7])
    if tested_header.get_sequence_version() != result_uniprot_lcr_headers[i][8]:
        print(tested_header.get_sequence_version(), "!=", result_uniprot_lcr_headers[i][8])
    if tested_header.get_lcr_start() != result_uniprot_lcr_headers[i][9]:
        print(tested_header.get_lcr_start(), "!=", result_uniprot_lcr_headers[i][9])
    if tested_header.get_lcr_end() != result_uniprot_lcr_headers[i][10]:
        print(tested_header.get_lcr_end(), "!=", result_uniprot_lcr_headers[i][10])


print('Pfam LCR headers...')
for i, header in enumerate(eg_pfam_lcr_headers):
    tested_header = st.PfamHeaderLCR(header)
    if tested_header.get_entry_name() != result_pfam_lcr_headers[i][0]:
        print(tested_header.get_entry_name(), "!=", result_pfam_lcr_headers[i][0])
    if tested_header.get_prot_start() != result_pfam_lcr_headers[i][1]:
        print(tested_header.get_prot_start(), "!=", result_pfam_lcr_headers[i][1])
    if tested_header.get_prot_end() != result_pfam_lcr_headers[i][2]:
        print(tested_header.get_prot_end(), "!=", result_pfam_lcr_headers[i][2])
    if tested_header.get_unique_identifer() != result_pfam_lcr_headers[i][3]:
        print(tested_header.get_unique_identifer(), "!=", result_pfam_lcr_headers[i][3])
    if tested_header.get_family_id() != result_pfam_lcr_headers[i][4]:
        print(tested_header.get_family_id(), "!=", result_pfam_lcr_headers[i][4])
    if tested_header.get_family_name() != result_pfam_lcr_headers[i][5]:
        print(tested_header.get_family_name(), "!=", result_pfam_lcr_headers[i][5])
    if tested_header.get_lcr_start() != result_pfam_lcr_headers[i][6]:
        print(tested_header.get_lcr_start(), "!=", result_pfam_lcr_headers[i][6])
    if tested_header.get_lcr_end() != result_pfam_lcr_headers[i][7]:
        print(tested_header.get_lcr_end(), "!=", result_pfam_lcr_headers[i][7])

print('Pfam in SP LCR headers...')
for i, header in enumerate(eg_pfam_in_sp_lcr_headers):
    tested_header = st.PfamInSPHeaderLCR(header)
    if tested_header.get_entry_name() != result_pfam_in_sp_lcr_headers[i][0]:
        print(tested_header.get_entry_name(), "!=", result_pfam_in_sp_lcr_headers[i][0])
    if tested_header.get_prot_start() != result_pfam_in_sp_lcr_headers[i][1]:
        print(tested_header.get_prot_start(), "!=", result_pfam_in_sp_lcr_headers[i][1])
    if tested_header.get_prot_end() != result_pfam_in_sp_lcr_headers[i][2]:
        print(tested_header.get_prot_end(), "!=", result_pfam_in_sp_lcr_headers[i][2])
    if tested_header.get_unique_identifer() != result_pfam_in_sp_lcr_headers[i][3]:
        print(tested_header.get_unique_identifer(), "!=", result_pfam_in_sp_lcr_headers[i][3])
    if tested_header.get_family_id() != result_pfam_in_sp_lcr_headers[i][4]:
        print(tested_header.get_family_id(), "!=", result_pfam_in_sp_lcr_headers[i][4])
    if tested_header.get_family_name() != result_pfam_in_sp_lcr_headers[i][5]:
        print(tested_header.get_family_name(), "!=", result_pfam_in_sp_lcr_headers[i][5])
    if tested_header.get_organism_name() != result_pfam_in_sp_lcr_headers[i][6]:
        print(tested_header.get_organism_name(), "!=", result_pfam_in_sp_lcr_headers[i][6])
    if tested_header.get_organism_identifier() != result_pfam_in_sp_lcr_headers[i][7]:
        print(tested_header.get_organism_identifier(), "!=", result_pfam_in_sp_lcr_headers[i][7])
    if tested_header.get_lcr_start() != result_pfam_in_sp_lcr_headers[i][8]:
        print(tested_header.get_lcr_start(), "!=", result_pfam_in_sp_lcr_headers[i][8])
    if tested_header.get_lcr_end() != result_pfam_in_sp_lcr_headers[i][9]:
        print(tested_header.get_lcr_end(), "!=", result_pfam_in_sp_lcr_headers[i][9])

print('Checking headers completed!')
print('^^ encountered errors above ^^')
