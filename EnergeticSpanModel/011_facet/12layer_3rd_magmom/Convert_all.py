from ase.io import read, write
import os

def convert_xsd_to_poscar():
  for file in os.listdir('.'):
    if file.lower().endswith('.xsd'):
      base_name = os.path.splitext(file)[0]
      poscar_name = base_name + '.poscar'

      try:
        atoms = read(file)
        write(poscar_name, atoms)
      except Exception as e:
          print(f"format transform failed : {file} ({e})")

def convert_poscar_to_xsd():
  for file in os.listdir('.'):
    if file.lower().endswith('.poscar'):
      base_name = os.path.splitext(file)[0]
      xsd_name = base_name + '.xsd'

      try:
        atoms = read(file)
        write(xsd_name, atoms)
      except Exception as e:
          print(f"format transform failed : {file} ({e})")


convert_xsd_to_poscar()
#convert_poscar_to_xsd()

