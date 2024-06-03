import pandas as pd
import pathlib
import vcf_reader


if __name__ == '__main__':
    results_folders = pathlib.Path('').rglob('*.vcf.gz')
    illumina_folder = pathlib.Path('ILLUMINA')
    nanopore_folder = pathlib.Path('NANOPORE')
    unique_files = sorted(set((file_path.name for file_path in results_folders)))
    output_folder = pathlib.Path('cartella_excel')

    for result in unique_files:
        vcf_path_nanopore = nanopore_folder / result
        vcf_path_illumina = illumina_folder / result

        print('='*100, result, '='*100, sep='\n')
        # print(vcf_path_illumina, vcf_path_nanopore, sep='\n'); continue
        df_nano = vcf_reader.dataframe(filename=vcf_path_nanopore, large=False)
        print(df_nano)
        df_illu = vcf_reader.dataframe(filename=vcf_path_illumina, large=False)
        print(df_illu)

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

        output_folder.mkdir(exist_ok=True)
        excel_path = output_folder / f'{result}_vcf.xlsx'
        df_differenze.to_excel(excel_path)
        print(f"File succesfully written: {excel_path}")
