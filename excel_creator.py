import pandas as pd
import pathlib
import vcf_handler


if __name__ == '__main__':
    results_folders = pathlib.Path('Output').iterdir()

    for result in results_folders:
        vcf_path_nanopore = (result / f'{result.name}_raw_and_isec_vcfs/0000.vcf.gz').as_posix()
        vcf_path_illumina = (result / f'{result.name}_raw_and_isec_vcfs/0001.vcf.gz').as_posix()

        # print('='*100, result, '='*100, sep='\n')
        df_nano = vcf_handler.dataframe(filename=vcf_path_nanopore, large=False)
        df_illu = vcf_handler.dataframe(filename=vcf_path_illumina, large=False)

        for i, df in enumerate([df_nano, df_illu]):
            colonne = df.columns.tolist()
            nome_colonna_tecnologia = 'TECNOLOGIA'
            colonne.insert(0, nome_colonna_tecnologia)
            if i == 0:
                df[nome_colonna_tecnologia] = 'NANOPORE'
                df_nano = df[colonne]
            else:
                df[nome_colonna_tecnologia] = 'ILLUMINA'
                df_illu = df[colonne]

        df_differenze = pd.concat([df_nano, df_illu])
        # print(df_differenze)

        excel_path = result / f'{result.name}_differenze_vcf.xlsx'
        df_differenze.to_excel(excel_path)
        print(f"File succesfully written: {excel_path}")
