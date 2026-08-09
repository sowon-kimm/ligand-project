from ase.io import read, write
atoms = read('CONTCAR')
write('ola-Hoxo_cont.xsd',atoms)
