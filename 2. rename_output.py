import pathlib, os

def rinomina(folder_path):
    gen = folder_path.rglob('README.txt')
    for result in gen:
        readme_row = 4
        isec_folder = result.parent
        with open(f'{isec_folder / "README.txt"}') as file:
            raw_content = file.readlines()

        for file in isec_folder.rglob('*.gz'):
            content = raw_content[readme_row]

            file_path, new_name, origin = content.split('\t')
            file_path = pathlib.Path(file_path)
            origin = origin.strip()
            origin = pathlib.Path(origin)

            if readme_row in [4, 5]:
                new_name = new_name.removeprefix('for ').replace(' ', '_')
                origin = f"{origin.parent.name}_{origin.name}"
                new_name = f"{new_name}_{origin}"
            else:
                new_name = new_name.removeprefix('for records from')
                new_name = new_name.removesuffix(' shared by both')
                new_name = pathlib.Path(new_name)
                new_name = f"{new_name.parent.name}_{new_name.name}"
                new_name = f"shared_records_in_notation_{new_name}"

            command = f"mv {file_path} {file_path.parent / new_name}"
            os.system(command)
            readme_row += 1

def pulisci(folder_path):
    gen = folder_path.rglob('README.txt')
    for result in gen:
        isec_folder = result.parent
        for index_file in isec_folder.rglob('*.gz.tbi'):
            os.system(f'rm "{index_file}"')

def estrai(folder_path):
    gen = folder_path.rglob('README.txt')
    for result in gen:
        isec_folder = result.parent
        os.system(f'rm "{result}"')
        for compressed_file in isec_folder.rglob('*.gz'):
            os.system(f'gunzip "{compressed_file}"')

if __name__ == '__main__':
    current_path = pathlib.Path('Output/').absolute()
    rinomina(current_path)
    pulisci(current_path)
    estrai(current_path)
