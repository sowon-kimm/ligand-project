from ase.io import read, write
atoms = read('CONTCAR')
write('olaoxo3_cont.xsd',atoms)
