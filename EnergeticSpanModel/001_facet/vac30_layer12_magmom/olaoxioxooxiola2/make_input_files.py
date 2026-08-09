from ase.calculators.vasp import Vasp
from ase.io import read, write
import os

os.environ['VASP_PP_PATH']="/home/shared/programs/vasp/vasp_pp"

input_path = "POSCAR"
atoms = read(input_path)

## fix atom ~14A
from ase.constraints import FixAtoms

mask = [atom.z < 16.5 for atom in atoms]
c = FixAtoms(mask=mask)
atoms.set_constraint(c)

## setting vasp calculation
vasp_calculator = Vasp(
    encut=520, xc='pbe', algo='fast', prec='Accurate',
    ismear=1, sigma=0.2, #electronic convergence
    ediff=1e-5, ispin=2,
    kpts=[4,4,1], gamma=True,
    ibrion=2, potim=0.5, ediffg=-0.05, isif=0, nsw=500, #geometric optimization
    nwrite=1, lcharg=False, lwave=False, lvtot=False, #output setting
    istart=0, lorbit=11,
    lreal='auto', npar=2,
    setups={'base': 'recommended'}, #paralleization
    ldau=True, 
    #ldau_luj={'Co': {'L': 2, 'U': 3.32, 'J': 0.0}}
    ldau_luj={'Mn': {'L': 2, 'U': 3.9, 'J': 0.0}}
)


## initialize
## set file name
#output_path = os.path.join(base_dir, folder_name)
#output_path=
#vasp_calculator.directory=output_path

vasp_calculator.write_input(atoms)
