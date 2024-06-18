import pandas as pd
import pathlib
import vcf_reader


CHOSEN_COLUMNS = ['INDICE_TECH', 'CAMPIONE', 'TECNOLOGIA', 'CHROM', 'POS', 'ID', 'REF',
                  'ALT', 'QUAL', 'FILTER', 'GT', 'AF', 'NVF', 'VF', 'GQ', 'DP', 'AD',
                  'DPD', 'DP4', 'PL', 'avsnp151', 'CLNALLELEID', 'CLNDN', 'CLNDISDB',
                  'CLNREVSTAT', 'CLNSIG']

if __name__ == '__main__':
    results_folders = pathlib.Path('./Output').iterdir()

    mega_df_list = []
    for result in results_folders:
        vcf_files = result.rglob('annotated*/result*.vcf')

        df_list = []
        for vcf_file in vcf_files:
            df = vcf_reader.dataframe(filename=vcf_file, large=False)
            df_columns = df.columns.tolist()
            technology_column_name = 'TECNOLOGIA'
            df_columns.insert(0, technology_column_name)
            if 'NANOPORE' in vcf_file.as_posix().upper():
                df[technology_column_name] = 'NANOPORE(4bases)'
            else:
                # df[technology_column_name] = 'ILLUMINA'
                if 'BRCA' in vcf_file.as_posix().upper() and not 'SOPHIA' in vcf_file.as_posix().upper():
                    df[technology_column_name] = 'ILLUMINA(Devyser)'
                else:
                    df[technology_column_name] = 'ILLUMINA(Sophia)'

            df = df[df_columns]
            df_list.append(df)
        
        try:
            differences_df = pd.concat(df_list)
        except ValueError:
            continue
        differences_df.index.name = 'INDICE_TECH'

        excel_columns = differences_df.columns.tolist()
        new_columns = [col for col in CHOSEN_COLUMNS if col in excel_columns]
        differences_df = differences_df[new_columns]
        if set([col for col in CHOSEN_COLUMNS if col in excel_columns]) != set(excel_columns).intersection(set(CHOSEN_COLUMNS)):
            print('colonne non filtrate correttamente')

        excel_path = result / f'{result.name}_differenze_vcf.xlsx'
        differences_df.to_excel(excel_path, index=True)
        print(f"File written: {excel_path}")
        differences_df['CAMPIONE'] = result.name
        new_columns.insert(0, 'CAMPIONE')
        differences_df.reset_index(drop=False, inplace=True)
        new_columns.insert(0, 'INDICE_TECH')
        differences_df = differences_df[new_columns]
        mega_df_list.append(differences_df)

    mega_df = pd.concat(mega_df_list)
    new_columns = [col for col in CHOSEN_COLUMNS if col in mega_df]

    mega_df = mega_df[new_columns]
    excel_path = pathlib.Path('Output/totale_differenze_vcf.xlsx')
    mega_df.to_excel(excel_path, index=False)
