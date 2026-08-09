from ase.io import read, write
atoms = read('CONTCAR')
write('1-2layer_cont.xsd',atoms)
