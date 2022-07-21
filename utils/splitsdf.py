#!/usr/bin/python

import sys,os
import argparse

 
parser = argparse.ArgumentParser(description='This script splits SDF file')
parser.add_argument('-i','--input', help='Input file name', required=True)
parser.add_argument('-n','--n_mol', help='Number of Molecules', required=True)
parser.add_argument('-o','--output', help='Output file name', required=True)
args = parser.parse_args()

infile = open(args.input)

n_mol = int(args.n_mol)

n = 0
n_file = 0

outfile = open('%s_0.sdf' %args.output, 'w')

for line in infile:
	outfile.write(line)
	if '$$$$' in line:
		n += 1
		if n % n_mol == 0:
			n_file += 1
			outfile = open('%s_%s.sdf' %(args.output, n_file), 'w')


sys.exit()

n_mol = int(args.n_mol)

mol_start_lines = []

with open(args.input) as infile:
	n = 0
	for line in infile:
		if 'V2000' in line:
			mol_start_lines.append(n-3)
		n += 1

infile = open(args.input)

print(mol_start_lines)

n_line = 0
i_mol = 0
n_file = 0
for line in infile:

	if i_mol % n_mol == 0:
		print('new file ', n_file)
		outfile = open('%s_%s.sdf' %(args.output, n_file), 'w')
		n_file += 1

	if n_line == mol_start_lines[i_mol]:
		i_mol += 1

	outfile.write(line)
	n_line += 1