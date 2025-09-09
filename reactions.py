s = ''

# hydrolysis
s += 'Hydrolysis\n'
X = ['COOH','Cl','NO3','OAc','OH']
for x in X:
    s += f"Mn({x})2 + H2O <-> Mn({x})(OH) + H{x}\n"

for x in X:
    s += f"Mn({x})(OH) + H2O <-> Mn(OH)2 + H{x}\n"

# ligand exchange
s += 'Ligand Exchange\n'
for x in X:
    for y in X:
        s += f"Mn({x})2 + H({y}) <-> Mn({x})({y}) + H{x}\n"

for x in X:
    for y in X:
        s += f"Mn({x})({y}) + H({y}) <-> Mn({y})2 + H{x}\n"

s += 'OAm reactions\n'

for x in X:
    s += f"OAm + H{x} <-> OAmH+ + {x}-\n"

s += 'Coordinate Covalent Bond\n'
for x in X:
    s += f"Mn({x})2 + 4H2O <-> Mn({x})2(OH2)4\n"

for x in X:
    s += f"Mn({x})2 + 4OAm <-> Mn({x})2(OH2)4\n"

for x in X:
    for y in X:
        s += f"Mn({x})({y}) + 4H2O <-> Mn({x})({y})(OH2)4\n"

for x in X:
    for y in X:
        s += f"Mn({x})({y}) + 4OAm <-> Mn({x})({y})(OAm)4\n"

s += 'Olation\n'
for x in X:
    s += f"*-({x}) + Mn-OH <->*-O-Mn + H({x})\n"

s += f"*=O- + Mn-OH <->*-O-Mn + OH-\n"
s += f"*-OH2+ + Mn-OH <->*-O-Mn + H3O+\n"

s += 'Oxolation\n'
for x in X:
    s += f"*-({x}) + Mn-OH2 <->*-OH-Mn + H({x})\n"

s += f"*=O- + Mn-OH2 <->*-OH-Mn + OH-\n"
s += f"*-OH2+ + Mn-OH2 <->*-OH-Mn + H3O+\n"

s += 'Never happens\n'
for x in X:
    s += f"*-({x}) + Mn=O- <->*-O-Mn + O-2\n"

s += f"*=O- + Mn=O- <->*-O-Mn + O-2\n"
s += f"*-OH2+ + Mn=O- <->*-O-Mn + O-2\n"


s += 'Capping\n'

for x in X:
    s += f"*-OH+H({x}) <-> *-({x}) + H2O\n"
    s += f"*-OH2+ + H({x}) <-> *-({x}) + H3O+\n"
    s += f"*=O-+H({x}) <-> *-({x}) + OH-\n"

for x in X:
    s += f"*-OAm+ + H({x}) <-> *-OAmH + {x}-\n"

with open("reactions.txt", "w", encoding="utf-8") as f:
    f.write(s)