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

def annotate_file(file_path, verbose):
    try:
        new_path = file_path.parent / 'temp'
        new_path.mkdir(exist_ok=True)
        new_file_path = new_path / file_path.name
        new_file_path = new_file_path.with_suffix('')  # Remove the .gz suffix and save as .vcf
        bed_file_path = Path('./bed_files/HCS_region_map.bed.gz').absolute()

        # command = f"zcat {file_path} | vcf-annotate -a {bed_file_path} -c CHROM,FROM,TO,INFO/ANN,INFO/Gene_Name -d key=INFO,ID=ANN,Number=1,Type=String,Description=Gene_Name -d key=INFO,ID=Gene_Name,Number=1,Type=String,Description=Gene_Name > {new_file_path}"
        # subprocess.run(command.split(' '), stdout=subprocess.DEVNULL)#, stderr=subprocess.DEVNULL)
        os.system(f"zcat {file_path} | \
                  vcf-annotate -a {bed_file_path} \
                  -c CHROM,FROM,TO,INFO/ANN,INFO/Gene_Name \
                  -d key=INFO,ID=ANN,Number=1,Type=String,Description='Gene Name' \
                  -d key=INFO,ID=Gene_Name,Number=1,Type=String,Description='Gene Name' \
                  > {new_file_path}")
        if verbose:
            print(f'Annotated: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error annotating {file_path}: {e}')

def filter_brca_file_genes(file_path, verbose):
    try:
        temp_file_path = file_path.with_suffix('').as_posix()
        temp_file_path = f"{temp_file_path}_temp.vcf"

        os.system(f'bcftools view -i \'INFO/Gene_Name="BRCA1" || INFO/Gene_Name="BRCA2"\' {file_path} > {temp_file_path}')

        os.system(f'cat {temp_file_path} > {file_path}')
        os.system(f'rm {temp_file_path}')
        if verbose:
            print(f'Filtered BRCA: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error filtering BRCA {file_path}: {e}')

def filter_file_genes(file_path, verbose):
    try:
        temp_file_path = file_path.with_suffix('').as_posix()
        temp_file_path = f"{temp_file_path}_temp.vcf"

        # bcftools view -i 'INFO/Gene_Name="ABRAXAS1" || INFO/Gene_Name="APC" || INFO/Gene_Name="ATM" || INFO/Gene_Name="BARD1" || INFO/Gene_Name="BRCA1" || INFO/Gene_Name="BRCA2" || INFO/Gene_Name="BRIP1" || INFO/Gene_Name="CDH1" || INFO/Gene_Name="CHEK2" || INFO/Gene_Name="EPCAM" || INFO/Gene_Name="MLH1" || INFO/Gene_Name="MRE11" || INFO/Gene_Name="MSH2" || INFO/Gene_Name="MSH6" || INFO/Gene_Name="MUTYH" || INFO/Gene_Name="NBN" || INFO/Gene_Name="PALB2" || INFO/Gene_Name="PIK3CA" || INFO/Gene_Name="PMS2" || INFO/Gene_Name="PMS2CL" || INFO/Gene_Name="PTEN" || INFO/Gene_Name="RAD50" || INFO/Gene_Name="RAD51C" || INFO/Gene_Name="RAD51D" || INFO/Gene_Name="STK11" || INFO/Gene_Name="TP53" || INFO/Gene_Name="XRCC2"' renamed_nanopore_annotated.vcf > renamed_nanopore_filtrato.vcf
        os.system(f'bcftools view -i \'INFO/Gene_Name="APC" || INFO/Gene_Name="ATM" || \
                  INFO/Gene_Name="BARD1" || INFO/Gene_Name="BRCA1" || INFO/Gene_Name="BRCA2" || INFO/Gene_Name="BRIP1" || \
                  INFO/Gene_Name="CDH1" || INFO/Gene_Name="CHEK2" || INFO/Gene_Name="EPCAM" || INFO/Gene_Name="MLH1" || \
                  INFO/Gene_Name="MSH2" || INFO/Gene_Name="MSH6" || INFO/Gene_Name="MUTYH" || INFO/Gene_Name="NBN" || \
                  INFO/Gene_Name="PALB2" || INFO/Gene_Name="PIK3CA" || INFO/Gene_Name="PMS2" || INFO/Gene_Name="PMS2CL" || \
                  INFO/Gene_Name="PTEN" || INFO/Gene_Name="RAD50" || INFO/Gene_Name="RAD51C" || INFO/Gene_Name="RAD51D" || \
                  INFO/Gene_Name="STK11" || INFO/Gene_Name="TP53" || INFO/Gene_Name="XRCC2"\' \
                  {file_path} > {temp_file_path}')

        os.system(f'cat {temp_file_path} > {file_path}')
        os.system(f'rm {temp_file_path}')
        if verbose:
            print(f'Filtered: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error filtering {file_path}: {e}')

def final_brca_filter(file_path, verbose):
    try:
        temp_file_path = file_path.as_posix().split('.')[0]
        temp_file_path = f"{temp_file_path}_temp.vcf"
        brca_bed_file_path = Path('./bed_files/BRCA_genes.bed').absolute()

        os.system(f'bcftools view -R {brca_bed_file_path} {file_path} > {temp_file_path}')

        os.system(f'cat {temp_file_path} > {file_path.with_suffix('')}')
        os.system(f'rm {temp_file_path}')
        if verbose:
            print(f'Filtered: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error filtering {file_path}: {e}')

def final_hc_filter(file_path, verbose):
    try:
        temp_file_path = file_path.with_suffix('').as_posix()
        temp_file_path = f"{temp_file_path}_temp.vcf"
        hc_bed_file_path = Path('./bed_files/NANOPORE_HC_BED_hg19.bed').absolute()

        os.system(f'bcftools view -R {hc_bed_file_path} {file_path} > {temp_file_path}')

        os.system(f'cat {temp_file_path} > {file_path.with_suffix('')}')
        os.system(f'rm {temp_file_path}')
        if verbose:
            print(f'Filtered: {file_path}')
    except subprocess.CalledProcessError as e:
        print(f'Error filtering {file_path}: {e}')

def get_matching_files(folder1, folder2, verbose):
    vcf_files = [file_path for file_path in folder1.rglob("*.vcf") if 'temp' in file_path.as_posix()]

    corresponding_paths = dict()
    no_matched_files = []
    for vcf_file in vcf_files:
        try:
            corresponding_paths[vcf_file] = next(folder2.rglob(vcf_file.name))
        except StopIteration:
            no_matched_files.append(vcf_file)

    if verbose:
        if len(no_matched_files) > 0:
            print('='*100)
            print('Missing matching file for the following files:')
            print('='*100)
            for file in no_matched_files:
                print(f"\t{file.as_posix()}")
        else:
            print('All files found a match!')

    return corresponding_paths, no_matched_files

def compare_files(comparison_tuples, output_dir, verbose):
    try:
        for file_tuple in comparison_tuples:
            file_comparison_name = file_tuple[0].name.split('.')[0]
            output_file_dir = output_dir / file_comparison_name
            output_file_dir.mkdir(exist_ok=True, parents=True)
            # bcftools stats renamed_nanopore_filtrato.vcf.gz illumina_filtrato.vcf.gz > results.vchk
            os.system(f'bcftools stats -c all {file_tuple[0]} {file_tuple[1]} > {output_file_dir / f"{output_file_dir.name}_comparison_stats.vchk"}')
            # vcf-compare renamed_nanopore_filtrato.vcf.gz illumina_filtrato.vcf.gz > compare.txt
            os.system(f'vcf-compare {file_tuple[0]} {file_tuple[1]} > {output_file_dir / f"{output_file_dir.name}_vcf_comparison.txt"}')
            # bcftools isec  -p isec_output -Oz illumina_filtrato.vcf.gz renamed_nanopore_filtrato.vcf.gz
            os.system(f'bcftools isec -c all -p {output_file_dir / f"{output_file_dir.name}_raw_and_isec_vcfs"} -Oz {file_tuple[0]} {file_tuple[1]}')
            if verbose:
                print(f'Comparison done for: {file_comparison_name}')
    except subprocess.CalledProcessError as e:
        print(f'Error comparing {file_tuple}: {e}')

def main(nanopore_folder_name, illumina_folder_name, 
         start_path = Path('.').absolute(), output_folder = Path('./Output').absolute(), 
         skip_confirmation = False, verbose = True):
    if not start_path.exists():
        raise ValueError('Invalid path provided')

    if not skip_confirmation:
        ask = input(f'Would you like to search in the folder {start_path}? (y/n): ')
        if ask.lower().strip() not in ['y', '']:
            exit(0)

    gz_file_generator = list_files_and_folders(start_path=start_path, extension='.gz')
    

    for file_path in gz_file_generator:
        recompress_file(file_path, verbose)
        index_file(file_path, verbose)
        str_upper_file_path = file_path.as_posix().upper()
        if ('BRCA' in str_upper_file_path or ('HC' in str_upper_file_path and nanopore_folder_name.upper() in str_upper_file_path)) and not ('ILLUMINA/BRCA_SOPHIA' in str_upper_file_path):
            normalize_file_chromosome_names(file_path, verbose)
            annotate_file(file_path, verbose)
        else:
            annotate_file(file_path, verbose)

    temp_vcf_files = list(file_path for file_path in list_files_and_folders(extension='.vcf') if 'temp' in file_path.as_posix())
    for file_path in temp_vcf_files:
        # filter_brca_file_genes(file_path, verbose)
        if 'BRCA' in file_path.as_posix().upper():
            filter_brca_file_genes(file_path, verbose)
            compress_file(file_path, verbose)
            file_path = Path(f"{file_path.as_posix()}.gz")
            index_file(file_path, verbose)
            final_brca_filter(file_path, verbose)
        else:
            filter_file_genes(file_path, verbose)
            compress_file(file_path, verbose)
            file_path = Path(f"{file_path.as_posix()}.gz")
            index_file(file_path, verbose)
            final_hc_filter(file_path, verbose)
    
    if verbose:
        print('='*50)
        print('Preprocessing Done!')
        print('='*50)

    # Getting matching names for each .vcf file
    nanopore_path = Path(start_path / nanopore_folder_name)
    illumina_path = Path(start_path / illumina_folder_name)
    matches_dict, missing_files = get_matching_files(nanopore_path, illumina_path, verbose)

    if len(missing_files) > 0:
        raise FileNotFoundError('Missing files! not all of the files have a match!')
    
    for file_path in temp_vcf_files:
        compress_file(file_path, verbose)
        file_path = f"{file_path.as_posix()}.gz"
        index_file(file_path, verbose)

    comparison_tuples = [(file1.with_name(file1.name + '.gz'), file2.with_name(file2.name + '.gz')) for file1, file2 in matches_dict.items()]
    compare_files(comparison_tuples, output_folder, verbose)


if __name__ == '__main__':
    main(nanopore_folder_name='NANOPORE', illumina_folder_name='ILLUMINA', 
         start_path=Path('.').absolute(), output_folder=Path('./Output').absolute(), 
         skip_confirmation=True, verbose=True)
