import pandas as pd
import pathlib
import vcf_reader


if __name__ == '__main__':
    results_folders = pathlib.Path('Output').iterdir()

    for result in results_folders:
        data_path = (result / f'{result.name}_raw_and_isec_vcfs')
        vcf_files = data_path.glob('*private*')

        df_list = []
        for vcf_file in vcf_files:
            df = vcf_reader.dataframe(filename=vcf_file, large=False)
            df_columns = df.columns.tolist()
            technology_column_name = 'TECNOLOGIA'
            df_columns.insert(0, technology_column_name)
            if vcf_file.as_posix().upper().endswith('NANOPORE.VCF'):
                df[technology_column_name] = 'NANOPORE'
            else:
                df[technology_column_name] = 'ILLUMINA'
            df = df[df_columns]
            df_list.append(df)
        
        differences_df = pd.concat(df_list)
        differences_df.index.name = 'INDICE_TECH'

        excel_path = result / f'{result.name}_differenze_vcf.xlsx'
        differences_df.to_excel(excel_path, index=True)
        print(f"File succesfully written: {excel_path}")
