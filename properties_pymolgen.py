#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 13:07:00 2021

@author: alexehaywood
"""

# Import Openeye Modules
import openeye.oechem as oe
from openeye import oemolprop as mp

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Are any H atoms attached to the atom?
def has_hydrogen(atom):
    if atom.GetImplicitHCount() > 0:
        return True
    for neigh in atom.GetAtoms():
        if neigh.IsHydrogen():

            return True

    return False

# Define Lipinski Acceptor
class IsLipinskiAcceptor(oe.OEUnaryAtomPred):
    def __call__(self, atom):
        if atom.GetAtomicNum() in [oe.OEElemNo_O, oe.OEElemNo_N]:

            return True

        return False

# Count the number of Hydrogen-Bond Acceptors 
def num_lipinsky_acceptors(mol):

    num_acc  =  oe.OECount(mol, IsLipinskiAcceptor())

    return num_acc

# Define Lipinsky Donor 
class IsLipinskiDonor(oe.OEUnaryAtomPred):
    def __call__(self, atom):
        if atom.GetAtomicNum() not in [oe.OEElemNo_O, oe.OEElemNo_N]:

            return False

        return has_hydrogen(atom)
    
# Count the number of Donors 
def num_lipinsky_donors(mol):

    num_don  =  oe.OECount(mol, IsLipinskiDonor())

    return num_don

# Calculate LogP and PSA
def oeMolProp(mol):

    props = [mp.OEGetXLogP(mol,atomxlogps = None), mp.OEGet2dPSA(mol,atomPSA = None)]

    return props

# Count the number of N atoms
def num_nitrogen_atom(mol_file):

    oe_in_file = oe.oemolistream()
    oe_in_file.open(mol_file)
    mols = oe_in_file.GetOEMols()
    num_nitrogen  =  [oe.OECount(mol, oe.OEIsNitrogen())  for mol in mols]

    return num_nitrogen[0]

# Count the number of O atoms
def num_oxygen_atom(mol_file):

    oe_in_file = oe.oemolistream()
    oe_in_file.open(mol_file)
    mols = oe_in_file.GetOEMols()
    num_oxygen  =  [oe.OECount(mol, oe.OEIsOxygen())  for mol in mols]

    return num_oxygen[0]

# Count the number of S atoms
def num_sulfur_atom(mol_file):

    oe_in_file = oe.oemolistream()
    oe_in_file.open(mol_file)
    mols = oe_in_file.GetOEMols()
    num_sulfur  =  [oe.OECount(mol, oe.OEIsSulfur())  for mol in mols]

    return num_sulfur[0]

# Count the number of C atoms
def num_carbon_atom(mol_file):

    oe_in_file = oe.oemolistream()
    oe_in_file.open(mol_file)
    mols = oe_in_file.GetOEMols()
    num_carbon  =  [oe.OECount(mol, oe.OEIsCarbon())  for mol in mols]

    return num_carbon[0]

# DEPRECATED: Count the number of atomatic systems 
def num_aromatic_ring_systems(mol_file):

    oe_in_file = oe.oemolistream()
    oe_in_file.open(mol_file)
    mols = oe_in_file.GetOEMols()
    nraromsystems  =  [oe.OEDetermineAromaticRingSystems(mol)  for mol in mols]
    
    return nraromsystems[0]

# Count the number of aromatic rings
def num_atomatic_rings(smi_mol):
     candidate = Chem.AddHs(Chem.MolFromSmiles(smi_mol))
     aroring = rdMolDescriptors.CalcNumAromaticRings(candidate)
     return aroring

# Count the number of rotatable bonds
def num_rot_bond(mol):

    n_rot_bonds  =  mp.OEGetRotatableBondCount(mol)

    return n_rot_bonds

# Count the number of chiral centers 
def num_chiral_centres(mol):
    
    n_chiral = oe.OECount(mol, oe.OEIsChiralAtom())
    
    return n_chiral

# Calculate Molecular Weight 
def molecular_weight(mol_file):

    oe_in_file = oe.oemolistream()
    oe_in_file.open(mol_file)
    mols = oe_in_file.GetOEMols()
    molw  =  [oe.OECalculateMolecularWeight(mol)  for mol in mols]

    return molw[0]

# Calculate sp3 fraction 
def sp3_fraction(smi_mol):
     candidate = Chem.AddHs(Chem.MolFromSmiles(smi_mol))
     csp3 = rdMolDescriptors.CalcFractionCSP3(candidate)
     return csp3