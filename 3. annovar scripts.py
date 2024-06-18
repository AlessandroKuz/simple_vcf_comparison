import pathlib
import os


if __name__ == '__main__':
    results_folders = pathlib.Path('./Output/').iterdir()

    table_annovar_path = pathlib.Path('/media/ngs_temp/annovar/table_annovar.pl')
    humandb_path = pathlib.Path('/media/ngs_temp/annovar/humandb')
    buildver = '--buildver hg19'
    protocol = '--protocol avsnp151,clinvar_20221231'
    operation = '--operation f,f'
    extra_params = '--remove --nastring . --polish --vcfinput'

    for result in results_folders:
        data_path = (result / f'{result.name}_raw_and_intersection_vcfs')
        vcf_files = data_path.glob('*private*')

        for vcf_file in vcf_files:
            output_name = vcf_file.parents[1].name
            output_dir_path = pathlib.Path(result / f'annotated_{output_name}_{"NANOPORE" if "NANOPORE" in vcf_file.as_posix() else "ILLUMINA"}')
            output_dir_path.mkdir(parents=False, exist_ok=True)
            output_param = f'--out {output_dir_path / "result"}'
            command = f'{table_annovar_path} {vcf_file} {humandb_path} {buildver} {output_param} {protocol} {operation} {extra_params}'
            os.system(command)
