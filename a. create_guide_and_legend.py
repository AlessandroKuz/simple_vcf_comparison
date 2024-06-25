from pathlib import Path


def write_guide(guide_path):
    guide = f"""{'='*100}\n{'\t'*6}Guida\n{'='*100}\n\n
Per ogni combinazione NANOPORE-ILLUMINA è presente una cartella con più file:\n
\t- FILE_vcf_comparison.txt  : (es: BRCA_118_21_vcf_comparison)
\t\t Contiene informazioni in formato testuale del diagramma di Venn dei record (quanti casi 
\t\t presenti solo in NANOPORE, quanti solo in ILLUMINA e quanti in comune)\n
\t- FILE_comparison_stats.txt: (es: BRCA_118_21_comparison_stats)
\t\t Contiene informazioni sulla quantità di record presenti per ogni file e la tipologia di mutazioni\n
\t- FILE_differenze_vcf.xlsx: (es: BRCA_118_21_differenze_vcf.xlsx)
\t\t Contiene i record delle differenze dei vcf con tecnologia e indice per ciascuno più le colonne 
\t\t possibilmente interessanti (descritte nella sottostante)\n
\t- FILE_raw_and_intersection_vcfs: (es: BRCA_118_21_raw_and_isec_vcfs)
\t\t Contiene i singoli vcf contenti le informazioni dei record presenti solo in NANOPORE, solo in ILLUMINA e quelli in comune\n


{'-'*100}\n{'\t'*6}Legenda\n{'-'*100}\n\n
TECNOLOGIA: Tecnologia impiegata nel sequenziamento
INDICE_TECH: Indice per i record di una determinata tecnologia
CHROM: -
POS: -
ID: -
REF: -
ALT: -
QUAL: -
FILTER: -
GT: Genotipo
AF: "Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall"
NVF: Probabilmente è AF Normalizzata (Normalized Variant Frequency) - presente solo nei SOPHIA, non sempre presente in tutti i record
VF: Variant Frequency stessa cosa di AF
GQ: "Genotype Quality"
DP: "Read depth at this position for the sample"
AD: "Allelic depth for the ref and alt alleles in the order listed"
DPD: "Read depth details: #ref plus strand, #ref minus strand, #alt plus strand, #alt minus strand"
DP4: "#ref plus strand, #ref minus strand, #alt plus strand, #alt minus strand"
PL: "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">

"""
    with open(guide_path, 'w') as file:
        file.write(guide)

if __name__ == '__main__':
    write_guide('Output/guida.txt')
