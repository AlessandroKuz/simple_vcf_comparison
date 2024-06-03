import re
import pprint


if __name__ == '__main__':
    genes = ['APC', 'ATM', 'BARD1', 'BRCA1', 'BRCA2', 'BRIP1', 'CDH1', 'CHEK2',
    'EPCAM', 'MLH1', 'MSH2', 'MSH6', 'MUTYH', 'NBN', 'PALB2', 'PIK3CA', 'PMS2',
    'PMS2CL', 'PTEN', 'RAD50', 'RAD51C', 'RAD51D', 'STK11', 'TP53', 'XRCC2']
    # genes = [f"{gene}\n" for gene in genes]
    
    with open('BED_RRMS_hg38.bed') as f:
        rows = f.readlines()

    # useful_rows = [row for row in rows if any(gene in row for gene in genes)]
    useful_rows = []
    for row in rows:
        for gene in genes:
            if row.endswith(f"\t{gene}\n"):
                row = row.lstrip('chr')
                useful_rows.append(row)
            # pattern = rf"\t{gene}\n$"
            # if re.findall(pattern, row):
            #     useful_rows.append(row)
    # if matches:
    #     print("Matches found:")
    #     for match in matches:
    #         print(match)
    print(len(useful_rows))
    useful_rows.sort(key=lambda x: int(x.split('\t')[0]))
    # # print(rows[0])
    # for char in rows[1000]:
    #     print([char], end='')
    # print(len([row for row in useful_rows if row.endswith('\n')]))

    pprint.pprint(useful_rows)
    with open('NANOPORE_HC_BED_hg19.bed', mode='w') as bed:
        bed.writelines(useful_rows)
