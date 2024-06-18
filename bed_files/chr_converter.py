import pandas as pd


def convert_chrom(in_bed_path, out_bed_path):
    with open(in_bed_path, 'r') as f:
        content = f.readlines()
    header = '\t'.join(['chr', 'start', 'stop', 'gene', '', '', '', '', '', '', '', ''])
    content.insert(0, header)
    with open(out_bed_path, 'w') as f:
        f.writelines(content)

    df = pd.read_table(in_bed_path, sep='\t', usecols=[0,1,2,3])
    df.iloc[:, 0] = df.iloc[:, 0].apply(lambda x: x.strip('chr'))
    df.to_csv(out_bed_path, sep='\t', index=False)

    with open(out_bed_path, 'r') as f:
        content = f.readlines()
    content.pop(0)
    with open(out_bed_path, 'w') as f:
        f.writelines(content)

if __name__ == '__main__':
    # convert_chrom('bed_files/raw/brca_hg38.bed', 'bed_files/brca_hg38_uscs.bed')
    # convert_chrom('bed_files/raw/hc_hg38.bed', 'bed_files/hc_hg38_uscs.bed')
    pass
