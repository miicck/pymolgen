import random

infile = open("/home/pczbf/Downloads/chembl_30_1m.sdf")

line_pos = []

n = 0
for line in infile:
	if 'RDKit' in line:
		line_pos.append(n-1)
	n += 1

print( line_pos[:10000])

infile.close()

infile = open("/home/pczbf/Downloads/chembl_30_1m.sdf")

random.seed(100)
random_pos = random.choice(line_pos)
print(random_pos)
lines = []

for i, line in enumerate(infile):
    if i >= random_pos:
        lines.append(line)
        if "$$$$" in line:
        	break

infile.close()

print(lines)
