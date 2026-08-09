from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.build import add_vacuum
import os

os.environ['VASP_PP_PATH']="/home/shared/programs/vasp/vasp_pp"

input_path = "CONTCAR"
atoms = read(input_path)

add_vacuum(atoms, vacuum = 10)

write('POSCAR', atoms, format='vasp')
