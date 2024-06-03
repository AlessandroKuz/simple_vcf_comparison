import pathlib, os

def rinomina(folder_path):
    gen = folder_path.rglob('README.txt')
    for result in gen:
        readme_row = 5
        isec_folder = result.parent
        # print(isec_folder)
        for file in isec_folder.rglob('*.gz'):
            # print(f"\t{file.relative_to(isec_folder)}")
            # print(file)
            # os.system(f'rm {file.as_posix() + ".tbi"}')
            new_name = os.popen(f"cat {isec_folder / 'README.txt'} | sed -n '{readme_row}p'")
            new_name = new_name.read().split('for records ')
            new_name = [element.strip('\t') for element in new_name]
            # file_name_to_rename = new_name[0]
            new_name = new_name [1:]
            new_name[0] = new_name[0].strip().strip('private to \t').strip('from ')
            # print("\t\t",(new_name))
            if readme_row in (5, 6):
                file_name = f'Records private to {file.parents[1].name}_{"ILLUMINA" if 'ILLUMINA' in new_name[0].upper() else "NANOPORE"}'
            else:
                file_name = f'Shared records from {file.parents[1].name}_{"ILLUMINA" if "ILLUMINA" in new_name[0].split(' ')[0] else "NANOPORE"}'
            # print(f'{file} --> {file_name}')
            # os.system(f'mv {file} \'{file.as_posix() + "/" + file_name + ".vcf.gz"}\'')
            os.system(f'mv {file} \'{file.parent.as_posix() + "/" + file_name + ".vcf.gz"}\'')
            # print('-'*150)
            readme_row += 1

def pulisci(folder_path):
    gen = folder_path.rglob('README.txt')
    for result in gen:
        isec_folder = result.parent
        # print(isec_folder)
        for index_file in isec_folder.rglob('*.gz.tbi'):
            os.system(f'rm "{index_file}"')

def estrai(folder_path):
    gen = folder_path.rglob('README.txt')
    for result in gen:
        isec_folder = result.parent
        # print(isec_folder)
        os.system(f'rm "{result}"')
        for compressed_file in isec_folder.rglob('*.gz'):
            os.system(f'gunzip "{compressed_file}"')

if __name__ == '__main__':
    current_path = pathlib.Path('Output/').absolute()
    rinomina(current_path)
    pulisci(current_path)
    estrai(current_path)
