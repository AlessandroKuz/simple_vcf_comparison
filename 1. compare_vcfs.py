from pathlib import Path
import subprocess, os


def list_files_and_folders(start_path='.', extension='.gz'):
    file_generator = (file_path for file_path in Path(start_path).rglob(f'*{extension}') 
                  if file_path.is_file() and 'bed' not in file_path.as_posix())
    return file_generator

def decompress_file(file_path, verbose):
    try:
        subprocess.run(['bgzip', '-d', file_path], check=True)
        if verbose:
            print(f'Decompressed: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error decompressing {file_path}: {e}')

def compress_file(file_path, verbose):
    try:
        subprocess.run(['bgzip', '-f', file_path], check=True)
        if verbose:
            print(f'Compressed: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error compressing {file_path}: {e}')

def recompress_file(file_path, verbose):
    decompress_file(file_path, verbose)
    compress_file(file_path.with_suffix(''), verbose)  # Remove the .gz suffix

def index_file(file_path, verbose):
    try:
        # bcftools index file_path
        subprocess.run(['bcftools', 'index', '-f', file_path], check=True)
        if verbose:
            print(f'Indexed File: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error indexing {file_path}: {e}')

def normalize_file_chromosome_names(file_path, verbose):
    try:
        # bcftools annotate --rename-chrs chromosome_nanopore.txt nanopore.vcf.gz -Oz -o renamed_nanopore.vcf.gz
        subprocess.run(['bcftools', 'annotate', '--rename-chrs', './maps/chromosome_mappings.txt', file_path.as_posix(), '-Oz', '-o', file_path.as_posix()], check=True)
        if verbose:
            print(f'Normalized CHROM: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error normalizing CHROM {file_path}: {e}')

def final_brca_filter(file_path, verbose):
    try:
        temp_file_path = file_path.with_suffix('').with_suffix('').as_posix()
        temp_file_path = f"{temp_file_path}_temp.vcf"
        # brca_bed_file_path = Path('./bed_files/brca_devyser_nanopore_hg19_chrless.bed')
        brca_bed_file_path = Path('./bed_files/Devyser_BRCA_hg19_chrless.bed')

        os.system(f'bcftools view -R {brca_bed_file_path} {file_path} > {temp_file_path}')

        os.system(f'cat {temp_file_path} > {file_path.with_suffix('')}')
        os.system(f'rm {temp_file_path} {file_path} {file_path.as_posix() + ".csi"}')
        compress_file(file_path.with_suffix(''), verbose)
        index_file(file_path, verbose)

        if verbose:
            print(f'Filtered: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error filtering {file_path}: {e}')

def final_hc_filter(file_path, verbose):
    try:
        temp_file_path = file_path.with_suffix('').with_suffix('').as_posix()
        temp_file_path = f"{temp_file_path}_temp.vcf"
        hc_bed_file_path = Path('./bed_files/hc_sophia_nanopore_hg19.bed')

        os.system(f'bcftools view -R {hc_bed_file_path} {file_path} > {temp_file_path}')

        os.system(f'cat {temp_file_path} > {file_path.with_suffix('')}')
        os.system(f'rm {temp_file_path} {file_path} {file_path.as_posix() + ".csi"}')
        compress_file(file_path.with_suffix(''), verbose)
        index_file(file_path, verbose)

        if verbose:
            print(f'Filtered: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error filtering {file_path}: {e}')

def get_matching_files(folder1, folder2):
    vcf_files = folder1.rglob("*.vcf.gz")

    folder2_type = folder2.name  # ILLUMINA or NANOPORE depending of what was passed as folder2
    corresponding_paths = dict()
    no_matched_files = []

    for vcf_file in vcf_files:
        equivalent_folder2_path = vcf_file.parents[1] / folder2_type / vcf_file.name
        if equivalent_folder2_path.exists():
            corresponding_paths[vcf_file] = equivalent_folder2_path
        else:
            no_matched_files.append(vcf_file)

    return corresponding_paths, no_matched_files

def compare_files(comparison_tuples, output_dir, verbose):
    try:
        for file_tuple in comparison_tuples:
            file_comparison_name = file_tuple[0].name.split('.')[0]
            output_file_dir = output_dir / file_comparison_name
            output_file_dir.mkdir(exist_ok=True, parents=True)
            # bcftools stats renamed_nanopore_filtrato.vcf.gz illumina_filtrato.vcf.gz > results.vchk
            os.system(f'bcftools stats -c all {file_tuple[0]} {file_tuple[1]} > {output_file_dir / f"{output_file_dir.name}_comparison_stats.txt"}')
            # vcf-compare renamed_nanopore_filtrato.vcf.gz illumina_filtrato.vcf.gz > compare.txt
            os.system(f'vcf-compare {file_tuple[0]} {file_tuple[1]} > {output_file_dir / f"{output_file_dir.name}_vcf_comparison.txt"}')
            # bcftools isec  -p isec_output -Oz illumina_filtrato.vcf.gz renamed_nanopore_filtrato.vcf.gz
            os.system(f'bcftools isec -c all -p {output_file_dir / f"{output_file_dir.name}_raw_and_intersection_vcfs"} -Oz {file_tuple[0]} {file_tuple[1]}')
            if verbose:
                print(f'Comparison done for: {file_comparison_name}')
    except subprocess.CalledProcessError as e:
        print(f'Error comparing {file_tuple}: {e}')

def main(nanopore_folder_name, illumina_folder_name, 
         start_path = Path('.').absolute(), output_folder = Path('./Output').absolute(), 
         verbose = True):
    if not start_path.exists():
        raise ValueError('Invalid path provided.')

    nanopore_path = Path(start_path / nanopore_folder_name)
    illumina_path = Path(start_path / illumina_folder_name)

    gz_file_generator_illumina = list_files_and_folders(start_path=illumina_path, extension='.gz')
    gz_file_generator_nanopore = list_files_and_folders(start_path=nanopore_path, extension='.gz')
    gz_files = list(gz_file_generator_illumina) + list(gz_file_generator_nanopore)
    gz_files.sort()

    if len(gz_files) == 0:
        raise FileNotFoundError('No vcf file was found.')

    for file_path in gz_files:
        recompress_file(file_path, verbose)
        index_file(file_path, verbose)
        str_upper_file_path = file_path.as_posix().upper()
        if ('BRCA' in str_upper_file_path or ('HC' in str_upper_file_path and nanopore_folder_name.upper() in str_upper_file_path)) and not ('ILLUMINA/BRCA_SOPHIA' in str_upper_file_path):
            normalize_file_chromosome_names(file_path, verbose)
        index_file(file_path, verbose)

    temp_vcf_files_illumina = list_files_and_folders(start_path=illumina_path, extension='.vcf.gz')
    temp_vcf_files_nanopore = list_files_and_folders(start_path=nanopore_path, extension='.vcf.gz')
    temp_vcf_files = list(temp_vcf_files_illumina) + list(temp_vcf_files_nanopore)
    temp_vcf_files.sort()

    for file_path in temp_vcf_files:
        if 'BRCA' in file_path.as_posix().upper() and not 'SOPHIA' in file_path.as_posix().upper():
            final_brca_filter(file_path, verbose)
        else:
            final_hc_filter(file_path, verbose)
    
    if verbose:
        print('='*50)
        print('Preprocessing Done!')
        print('='*50)

    # Getting matching names for each .vcf file
    matches_dict, missing_files = get_matching_files(illumina_path, nanopore_path)

    if len(matches_dict) == 0: 
        raise FileExistsError('No file to match!')
    if len(missing_files) > 0:
        raise FileNotFoundError('Missing files! Not all of the files have a match!')

    comparison_tuples = list(matches_dict.items())
    comparison_tuples.sort(key=lambda x: x[0])
    compare_files(comparison_tuples, output_folder, verbose)


if __name__ == '__main__':
    main(nanopore_folder_name='NANOPORE', illumina_folder_name='ILLUMINA', 
         start_path=Path('.').absolute(), output_folder=Path('./Output').absolute(), 
         verbose=True)
