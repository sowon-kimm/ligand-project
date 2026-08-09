import os
import shutil

source_dir = "./xsd_files/"

for filename in os.listdir(source_dir):
    if filename.endswith('.poscar'):
        folder_name = filename[:-7]
        target_dir = os.path.join('.', folder_name)       

        os.makedirs(target_dir, exist_ok=True)
        
        src = os.path.join(source_dir, filename)
        dst = os.path.join(target_dir, 'POSCAR')
	
#        file1 = 'make_input_files.py'
#        file2 = 'submit_vasp.sh'

        shutil.copy(src, dst)
#        shutil.copy(file1, os.path.join(target_dir, file1))
#        shutil.copy(file2, os.path.join(target_dir, file2))
        print(f"폴더 생성: {target_dir}")
