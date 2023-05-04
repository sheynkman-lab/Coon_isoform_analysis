#!/usr/bin/python3
import pandas as pd
from collections import defaultdict

still_testing = False

"""
Working with data from:
Global detection of human variants and isoforms by deep proteome sequencing.
Sinitcyn P, Richards AL, Weatheritt RJ, Brademan DR, Marx H, Shishkova E, Meyer JG, Hebert AS, Westphall MS, Blencowe BJ, Cox J, Coon JJ
Nat Biotechnol. 2023 Mar 23
DOI: 10.1038/s41587-023-01714-x, PMID: 36959352

Supplemental table (41587_2023_1714_MOESM3_ESM.xlsx) from: https://www.nature.com/articles/s41587-023-01714-x#Sec24
Gencode data (gencode.v43c_translations.fa): https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz
"""

coon = pd.read_excel('41587_2023_1714_MOESM3_ESM.xlsx', sheet_name ='Table S3', converters={'GeneName':str,})

# To see this, sort by Gene_name.  A number of entries will float to the top where the gene name gets converted to a date. 
coon_genes = coon['GeneName'].sort_values().to_list()
for gene in coon_genes[0:200]:
    if '2020' in gene:
        coon_genes.remove(gene)
print (f"coon_genes is {len(coon_genes)} records long") #print (coon_genes)

"""
## Build a dict of peptides 
peptides[gene] = [ list of peptide strings ] 
In the imported dataframe, if there are peptides the dtype is str, otherwise NaN
"""

all_peptides={}
columns =  ['GeneName', 'Proteomics|aspn|Peptides|Path0', 'Proteomics|aspn|Peptides|Path1', 
                             'Proteomics|chymo|Peptides|Path0', 'Proteomics|chymo|Peptides|Path1', 
                             'Proteomics|gluc|Peptides|Path0', 'Proteomics|gluc|Peptides|Path1', 
                             'Proteomics|lysc|Peptides|Path0', 'Proteomics|lysc|Peptides|Path1', 
                             'Proteomics|lysn|Peptides|Path0', 'Proteomics|lysn|Peptides|Path1', 
                             'Proteomics|trypsin|Peptides|Path0', 'Proteomics|trypsin|Peptides|Path1' ]
coon_red = coon.loc[:,columns]  # coon_red short for coon_reduced
for i,row in coon_red.iterrows():
    if row[0] not in all_peptides:
        all_peptides[row[0]] = []
    for col in range(1,13):
        if type(row[col]) == str :
            peps=row[col].split(';')
            all_peptides[row[0]].extend(peps)
print (f'{len(all_peptides)} genes are represented in all_peptides')

"""
    ##  Build data structure for genes
    genes[gene] = { isoform_name : sequences }
"""
''' Builds a data structure for genes. genes is a dict with gene names as keys 
    The 'values' of the genes[gene] is a dict with {isoform_name:sequence}
'''
gene_file='gencode.v43.pc_translations.fa'
genes ={}
with open(gene_file, "r") as f:
    lines = f.read().splitlines()
for line in lines:
    if line.startswith(">"):
        parts = line.split('|')
        gene = parts[6]
        isoform = parts[1]
        if gene not in genes:
            genes[gene] = {}# defaultdict()
        genes[gene][isoform] = ""
    else:
        genes[gene][isoform] += line


"""
# Big loop
The function gene_matches compares all the peptides against all of the translated protein isoforms in Gencode 43.  
If a peptide only matches 1 isoform, its ENST is added to 1) the set iso_matches and 2) the dict single_iso_matches_dict as {peptide:ENST}

#### Output looks like this if still_testing = True
  
1 peptides found for gene A1CF
	EIYMNVPVGAAGVR	matches ENST00000373993.6
	EIYMNVPVGAAGVR	matches ENST00000395489.7
		Peptide EIYMNVPVGAAGVR matched 2 isoforms. {'ENST00000395489.7', 'ENST00000373993.6'}
2 peptides found for gene A2ML1
	HLHCISFLVPPPAGGTEEVATIR	matches ENST00000299698.12
	YSMVELQDPNSNR	matches ENST00000299698.12
4 peptides found for gene AAAS
	KFAVALLDDSVRVYNASSTIVPSL	matches ENST00000548931.6
	KFAVALLDDSVRVYNASSTIVPSL	matches ENST00000209873.9
	KFAVALLDDSVRVYNASSTIVPSL	matches ENST00000550286.5
	SEDLIAEFAQVTNWSSCCLR	matches ENST00000209873.9
	SEDLIAEFAQVTNWSSCCLR	matches ENST00000550286.5
	FAVALLDDSVRVYNASSTIVPSLK	matches ENST00000548931.6
	FAVALLDDSVRVYNASSTIVPSLK	matches ENST00000209873.9
	FAVALLDDSVRVYNASSTIVPSLK	matches ENST00000550286.5
	VYNASSTIVPSLK	matches ENST00000548931.6
	VYNASSTIVPSLK	matches ENST00000209873.9
	VYNASSTIVPSLK	matches ENST00000550286.5
		Peptide KFAVALLDDSVRVYNASSTIVPSL matched 3 isoforms. {'ENST00000209873.9', 'ENST00000550286.5', 'ENST00000548931.6'}
		Peptide SEDLIAEFAQVTNWSSCCLR matched 2 isoforms. {'ENST00000209873.9', 'ENST00000550286.5'}
		Peptide FAVALLDDSVRVYNASSTIVPSLK matched 3 isoforms. {'ENST00000209873.9', 'ENST00000550286.5', 'ENST00000548931.6'}
		Peptide VYNASSTIVPSLK matched 3 isoforms. {'ENST00000209873.9', 'ENST00000550286.5', 'ENST00000548931.6'}
"""

def gene_matches(gene, still_testing=False):
    peptides =  set(all_peptides[gene])
    if still_testing:    
        print (f'{len(peptides)} peptides found for gene {gene}')
    try:
        proteins = genes[gene]
    except:
        return (0,0)
    #print (proteins)
    #Create a set of protein matches for each peptide
    protein_matches = defaultdict(set)
    for peptide in peptides:
        for iso,seq in proteins.items():
            if peptide in seq:
                if still_testing:
                    print (f'\t{peptide}\tmatches {iso}')
                protein_matches[peptide].add(iso)

    # Identify peptides that only match one protein
    iso_matches = set()
    single_iso_matches_dict = {}
    for peptide, iso_set in protein_matches.items():
        if len(iso_set) == 1:
            iso = iso_set.pop()
            single_iso_matches_dict[peptide] = iso
            iso_matches.add(iso)
        else:
            if still_testing:
                print (f'\t\tPeptide {peptide} matched {len(iso_set)} isoforms. {iso_set}')
    #print("Single protein matches dict:", single_iso_matches_dict)
    return (iso_matches, single_iso_matches_dict)


all_iso_matches = {}
all_single_iso_matches_dict = {}
if still_testing:
    to_check = coon_genes[0:5]
else:
    to_check = coon_genes
for gene in to_check:
#for gene in ['SDF4']":
    iso_matches, single_iso_matches_dict = gene_matches(gene, still_testing= still_testing)
    if iso_matches:
        all_iso_matches[gene] = iso_matches
        all_single_iso_matches_dict[gene] = single_iso_matches_dict
if still_testing:
    print (all_single_iso_matches_dict)
    print (all_iso_matches)

multiple_isoforms = []
for gene, isos in all_iso_matches.items():
    if len(isos) > 1:
        multiple_isoforms.append(gene)
        #print (f'Gene {gene} has more than one isoform!')
print (f'{len (multiple_isoforms)} genes have multiple isoforms that have unique peptides')        
print (multiple_isoforms)

