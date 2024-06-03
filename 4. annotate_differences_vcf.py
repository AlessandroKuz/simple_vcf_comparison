import pandas as pd
import hgvs
import hgvs.parser
import vcf_reader
import pathlib
import requests
import json
import pprint


def get_variant_info_grch37(hgvsg_not: str):
    url = f"https://grch37.rest.ensembl.org/vep/human/hgvs/{hgvsg_not}?content-type=application/json"
    response = requests.get(url, headers={"Content-Type": "application/json"})

    if response.status_code == 200:
        return json.loads(response.content.decode("utf-8"))
    else:
        print(f"Error: {response.status_code}")
        print(f"Response content: {response.content}")
        # if len(ref) > len(alt): del?
        # elif len(alt) > len(ref): ins, dup?
        # else MNPs?
        # for variazione_possibile in variazioni_possibili:
            # prova = ...
            # if prova:
                # return risultato

        return None

# def get_annotations(hgvsg_not: str) -> list:
#     """
#     Get the HGVS from a VCF file
#     :param hgvsg_not:
#     :return:
#     """
#     variants = []
#     callset = al.vcf_to_dataframe(vcf_path, fields='*')

#     def get_variant_info_grch37(hgvs: str):
#         url = f"https://grch37.rest.ensembl.org/vep/human/hgvs/{hgvs}?content-type=application/json"
#         response = requests.get(url, headers={"Content-Type": "application/json"})

#         if response.status_code == 200:
#             return json.loads(response.content.decode("utf-8"))
#         else:
#             print(f"Error: {response.status_code}")
#             print(f"Response content: {response.content}")
#             return None

#     for index, row in callset.iterrows():
#         variants.append(get_variant_info_grch37(row['HGVSG']))

#     return variants

if __name__ == '__main__':
    # vcf_path = pathlib.Path('ILLUMINA/BRCA_327_22.vcf.gz')
    # df = vcf_reader.dataframe(vcf_path, large=False)
    # print(df)
    parser_hgvs = hgvs.parser.Parser()
    # print(hgvs.posedit.PosEdit(df.loc[0, 'HGVSG']))
    # print(parser_hgvs(grch37))

    all_vcf_files = pathlib.Path('Output').rglob('*.vcf.gz')
    variant_annotations = dict()

    for vcf_file_path in all_vcf_files:
        df = vcf_reader.dataframe(vcf_file_path, large=False)
        try:
            df.loc[:, 'HGVSG']
            has_hgvsg_notation = True
        except KeyError:
            has_hgvsg_notation = False
            print(f"VCF {vcf_file_path} has no HGVSG notion, so not all entries might be correct.")

        for index in range(len(df)):
            if not has_hgvsg_notation:
                chrom = df.loc[index, 'CHROM']
                try:
                    int(chrom[0])
                    chrom = f'chr{chrom}'
                except ValueError:
                    pass

                pos = int(df.loc[index, 'POS'])
                ref = df.loc[index, 'REF']
                alt = df.loc[index, 'ALT']

                posedit = hgvs.posedit.PosEdit(
                        hgvs.location.Interval(hgvs.location.SimplePosition(pos), hgvs.location.SimplePosition(pos + len(alt) - 1)),
                        hgvs.edit.NARefAlt(ref=ref, alt=alt))
                simple_posedit = hgvs.posedit.PosEdit(hgvs.location.SimplePosition(pos), hgvs.edit.NARefAlt(ref, alt))
                
                print(str(posedit) == str(simple_posedit))
                var_g = str(chrom) + ':g.' + str(posedit)
            else:
                try:
                    hgvsg = df.loc[index, 'HGVSG']
                    int(hgvsg[0])
                    hgvsg = f'chr{hgvsg}'
                except ValueError:
                    pass
                var_g = df.loc[index, 'HGVSG']

            skippa = True if '>' in var_g else False
            if skippa:
                continue

            existing_value = variant_annotations.get(var_g)  # check for existing value
            if existing_value:
                continue  # avoids unnecessary calls for already existing variants
            var_annotation = get_variant_info_grch37(var_g)
            if var_annotation:
                var_annotation = var_annotation[0]
                pprint.pprint(var_annotation)
            variant_annotations[var_g] = var_annotation

            if has_hgvsg_notation and hgvsg != var_g:
                print(posedit)
                print(parser_hgvs.parse_hgvs_variant(var_g))
                print(parser_hgvs.parse_hgvs_variant(hgvsg))

    # crea qua qualcosa che scrive le info in un file json
    json.dump(variant_annotations, 'annotations.json')

# TODO: 
# [ ] Crea un colonna HGVSG per ogni file che poi fungerà da chiave di unione:
# [ ] Implementare logica per la conversione in del, ins, dup, delins, ...
# [X] Prima di fare le chiamate crea tutti gli HGVSG e controlla se variabile già esiste
# [X] Salva risultati in un file
# [ ] se vuoi ricontrolla le non annotate
# [ ] uplifting da grch37 a grch38 e utilizzo nuovo bed

'''
'47641559delinsT' --> '47641559insT'
'chr2:g.47641559delinsT' --> 'chr2:g.47641559delA'

'''

def vcf_to_hgvs(chrom, pos, ref, alt):
    # Construct the HGVS variant object
    if len(ref) == 1 and len(alt) == 1:
        var_hgvs = f"{chrom}:g.{pos}{ref}>{alt}"
    elif len(ref) > len(alt):
        var_hgvs = f"{chrom}:g.{pos}_{pos+len(ref)-1}del{ref[1:]}"
    elif len(ref) < len(alt):
        var_hgvs = f"{chrom}:g.{pos}_{pos+1}ins{alt[1:]}"
    else:
        var_hgvs = f"{chrom}:g.{pos}_{pos+len(ref)-1}delins{alt}"

"""
Abbiamo varie tipologie di mutazioni:
    1. SNP (Single Nucleotide Polymorphism):
        Description: A single nucleotide change.
        Identification: Both REF and ALT are single nucleotides of the same length (one nucleotide each).
        Example: 
            REF     ALT
            A       G
    2. MNP (Multi-Nucleotide Polymorphism):
        Description: A change involving multiple adjacent nucleotides.
        Identification: REF and ALT have the same length and both are longer than one nucleotide.
        Example: 
            REF = AGT, ALT = TCA.
    3. Insertion:
        Description: One or more nucleotides are inserted in the sequence.
        Identification: ALT is longer than REF (i.e., ALT sequence has extra nucleotides).
        Example: 
            REF = A, ALT = ATG (insertion of TG).
    4. Deletion:
        Description: One or more nucleotides are deleted from the sequence.
        Identification: REF is longer than ALT (i.e., REF sequence has extra nucleotides).
        Example: 
            REF = ATG, ALT = A (deletion of TG).
    5. Indel (Insertion/Deletion):
        Description: A mutation that involves both an insertion and a deletion, changing the sequence length.
        Identification: REF and ALT are of different lengths and involve both insertion and deletion components.
        Example: 
            REF = AGT, ALT = A (deletion of GT and possible insertion).
    6. Complex (or Block Substitution):
        Description: Multiple changes within a short region that cannot be classified as simple SNPs or indels.
        Identification: Both REF and ALT are longer than one nucleotide and differ in a non-uniform way.
        Example: 
            REF = AGT, ALT = TCG.
"""


"""
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import hgvs.variantmapper
import biocommons.seqrepo

# Initialize hgvs data provider and mapper
hdp = hgvs.dataproviders.uta.connect()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")

# Create an instance of the HGVS parser
hp = hgvs.parser.Parser()

# Function to convert VCF to HGVS g. notation
def vcf_to_hgvs(chrom, pos, ref, alt):
    # Construct the HGVS variant object
    if len(ref) == 1 and len(alt) == 1:
        var_hgvs = f"{chrom}:g.{pos}{ref}>{alt}"
    elif len(ref) > len(alt):
        var_hgvs = f"{chrom}:g.{pos}_{pos+len(ref)-1}del{ref[1:]}"
    elif len(ref) < len(alt):
        var_hgvs = f"{chrom}:g.{pos}_{pos+1}ins{alt[1:]}"
    else:
        var_hgvs = f"{chrom}:g.{pos}_{pos+len(ref)-1}delins{alt}"

    # Parse the HGVS string to create a Variant object
    var = hp.parse_hgvs_variant(var_hgvs)

    # Normalize and map to genomic coordinates
    var_g = am.c_to_g(var)

    return str(var_g)

# Example usage
chrom = "1"
pos = 123456
ref = "A"
alt = "G"
hgvs_notation = vcf_to_hgvs(chrom, pos, ref, alt)
print(hgvs_notation)  # Output: 1:g.123456A>G
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import hgvs.variantmapper
import biocommons.seqrepo

# Initialize hgvs data provider and mapper
hdp = hgvs.dataproviders.uta.connect()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")

# Create an instance of the HGVS parser
hp = hgvs.parser.Parser()

# Function to convert VCF to HGVS g. notation
def vcf_to_hgvs(chrom, pos, ref, alt):
    # Construct the HGVS variant object
    if len(ref) == 1 and len(alt) == 1:
        var_hgvs = f"{chrom}:g.{pos}{ref}>{alt}"
    elif len(ref) > len(alt):
        var_hgvs = f"{chrom}:g.{pos}_{pos+len(ref)-1}del{ref[1:]}"
    elif len(ref) < len(alt):
        var_hgvs = f"{chrom}:g.{pos}_{pos+1}ins{alt[1:]}"
    else:
        var_hgvs = f"{chrom}:g.{pos}_{pos+len(ref)-1}delins{alt}"

    # Parse the HGVS string to create a Variant object
    var = hp.parse_hgvs_variant(var_hgvs)

    # Normalize and map to genomic coordinates
    var_g = am.c_to_g(var)

    return str(var_g)

# Example usage
chrom = "1"
pos = 123456
ref = "A"
alt = "G"
hgvs_notation = vcf_to_hgvs(chrom, pos, ref, alt)
print(hgvs_notation)  # Output: 1:g.123456A>G
"""