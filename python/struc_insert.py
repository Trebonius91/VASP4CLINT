#!/usr/bin/env python3.6
import sys 
import os 
import numpy as np
from scipy.spatial import distance_matrix

print('''
This script inserts a given structure A (e.g., an intermetallic 
phase crystal) into another structure B (e.g., a SCALMS droplet).
For this, all atoms of B in the region where A is placed will be 
removed.
The unitcell of B should be much larger than A in order to contain it.
The dimensions of B will be taken as dimensions of the result.
Strukture A must be stored in POSCAR_insert, structure B must be stored
in POSCAR_solvent. Both must be given in cartesian coordinates!
In addition, the desired center of mass position of structure A must
be given in the file insert_com.dat (x,y,z in one line, in Angstrom)

The result will be written in POSCAR.
''')

#
#    Global, predefined parameters
#
#    the distance, below which atoms of the solvent are removed (in Angstrom)

dist_remove = 2.8


# Dictionary of all elements matched with their atomic masses.
elements_dict = {'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'Be' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,\
                 'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,\
                 'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,\
                 'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,\
                 'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,\
                 'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,\
                 'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,\
                 'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,\
                 'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,\
                 'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,\
                 'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,\
                 'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,\
                 'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'GD' : 157.25,\
                 'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,\
                 'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,\
                 'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,\
                 'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,\
                 'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,\
                 'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,\
                 'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,\
                 'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247, 'Bk' : 247,\
                 'Ct' : 251, 'Es' : 252, 'Fm' : 257, 'Md' : 258, 'No' : 259,\
                 'Lr' : 262, 'Rf' : 261, 'Db' : 262, 'Sg' : 266, 'Bh' : 264,\
                 'Hs' : 269, 'Mt' : 268, 'Ds' : 271, 'Rg' : 272, 'Cn' : 285,\
                 'Nh' : 284, 'Fl' : 289, 'Mc' : 288, 'Lv' : 292, 'Ts' : 294,\
                 'Og' : 294}


#
#    Read in the coordinates of structure A (inserted crystal or structure)
#
insert_name="POSCAR_insert"

insert_in = open(insert_name,"r")

with insert_in as infile:
   line = infile.readline()
   line = infile.readline()
   line = infile.readline()
   line = infile.readline()
   line = infile.readline()
   line = infile.readline()
#   The elements line
   line = line.rstrip("\n")
   elements_insert = line.rstrip().split()
   nelem_insert = len(elements_insert)
   line = infile.readline().rstrip("\n")
   line_split = line.rstrip().split()
   elem_num_insert=[]
   natoms_insert=0
   names_insert=[]
   for i in range(nelem_insert):
      elem_num_insert.append(int(line_split[i]))
      # total number of atoms in the surface
      natoms_insert=natoms_insert+elem_num_insert[i]
      for j in range(elem_num_insert[i]):
          names_insert.append(elements_insert[i])
  
   natoms_surf=int(natoms_insert)

   xyz_insert = np.zeros((natoms_insert,3))
   mass_insert = np.zeros(natoms_insert)
   line = infile.readline()
   line = line.rstrip("\n")
   print(line)
   if line == "Direct" or line == "direct":
      print(" Please give POSCAR_insert in cartesian coordinates!")
      sys.exit()
 
   for i in range(natoms_insert):
      line = infile.readline().rstrip("\n")
      xyz_read = line.rstrip().split()[0:3]
      for j in range(3):
         xyz_insert[i][j]=float(xyz_read[j])
#    Determine insertion center of mass


   com=np.zeros(3)
   mass=0.0
   for i in range(natoms_insert):
      for j in range(3):
         com[j]=com[j]+elements_dict[names_insert[i]]*xyz_insert[i][j]
      mass=mass+elements_dict[names_insert[i]]

   for i in range(3):
      com[i]=com[i]/mass

 
#
#    Read in the coordinates of structure B (solvent)
#

solvent_name="POSCAR_solvent"

solvent_in = open(solvent_name,"r")

with solvent_in as infile:
   line = infile.readline()
   line = infile.readline()
   line = infile.readline().rstrip("\n")
   unit_line1 = line
   line = infile.readline().rstrip("\n")
   unit_line2 = line
   line = infile.readline().rstrip("\n")
   unit_line3 = line
   line = infile.readline()

#   The elements line
   line = line.rstrip("\n")
   elements_solvent = line.rstrip().split()
   nelem_solvent = len(elements_solvent)
   line = infile.readline().rstrip("\n")
   line_split = line.rstrip().split()
   elem_num_solvent=[]
   natoms_solvent=0
   names_solvent=[]
   for i in range(nelem_solvent):
      elem_num_solvent.append(int(line_split[i]))
      # total number of atoms in the surface
      natoms_solvent=natoms_solvent+elem_num_solvent[i]
      for j in range(elem_num_solvent[i]):
          names_solvent.append(elements_solvent[i])

   natoms_solvent=int(natoms_solvent)

   xyz_solvent = np.zeros((natoms_solvent,3))
   line = infile.readline()
   line = line.rstrip("\n")

   if line == "Direct" or line == "direct":
      print(" Please give POSCAR_solvent in cartesian coordinates!")
      sys.exit()

   for i in range(natoms_solvent):
      line = infile.readline().rstrip("\n")
      xyz_read = line.rstrip().split()[0:3]
      for j in range(3):
         xyz_solvent[i][j]=float(xyz_read[j])

#
#    Read in the desired center of mass position of the inserted structure
#
com_name="insert_com.dat"

com_in = open(com_name,"r")

com_insert = np.zeros(3)

with com_in as infile:
   line = infile.readline()
   com_read = line.rstrip().split()
   for i in range (3):
      com_insert[i] = float(com_read[i])
      
   
#    Shift inserted structure to desired COM
com_shift=com_insert-com

for i in range(natoms_insert):
   for j in range(3):
      xyz_insert[i][j]=xyz_insert[i][j]+com_shift[j]


print(com_shift)

#
#     Compare both structures, generate a boolean mask of deleted atoms in B
#

dist_mat = distance_matrix(xyz_solvent,xyz_insert,p=2)

dist_mins = dist_mat.min(axis=1)

min_mask=dist_mins < dist_remove 

#
#    Now print combined structure with insertion to POSCAR: All atoms of 
#    structure A and atoms of structure B that are far enough away from A
#    A will be printed after B, with two separated element sections
#

original_stdout=sys.stdout
with open("POSCAR","w") as f:
   sys.stdout = f

   print("Insertion into solvent system, made with scalms_insert.py")
   print(1.0)
   print(unit_line1)
   print(unit_line2)
   print(unit_line3)
   
#
#    Determine number of element atoms in solvent
#
   xyz_solv_print=np.zeros((natoms_solvent,3))
   elem_num_solv_new=[0]*nelem_solvent
   natoms_solv_new=0    
   index=0
   index2=0 
#
#    Filter the deleted atoms of the solvent
# 
   for i in range(nelem_solvent):
      # total number of atoms in the surface
      for j in range(elem_num_solvent[i]):
          if not min_mask[index]:
             xyz_solv_print[index2][0]=xyz_solvent[index][0]
             xyz_solv_print[index2][1]=xyz_solvent[index][1]
             xyz_solv_print[index2][2]=xyz_solvent[index][2]
             index2=index2+1
             elem_num_solv_new[i]=elem_num_solv_new[i]+1
             natoms_solv_new=natoms_solv_new+1
          index=index+1
#
#    Print the element names section
#
   for i in range(nelem_solvent):
      if elem_num_solv_new[i] > 0:
         print(elements_solvent[i], end=" ")        
   for i in range(nelem_insert):
      print(elements_insert[i], end=" ")
   print(" ")

#
#    Print the element numbers section
#
  
   for i in range(nelem_solvent):
      if elem_num_solv_new[i] > 0:
         print(elem_num_solv_new[i], end=" ")
   for i in range(nelem_insert):
      print(elem_num_insert[i], end=" ")
   print(" ")
 
   print("Cartesian")
   
   for i in range(natoms_solv_new):
      print(str(xyz_solv_print[i][0]) + " " + str(xyz_solv_print[i][1]) + " " + 
           str(xyz_solv_print[i][2])) 

   for i in range(natoms_insert):
      print(str(xyz_insert[i][0]) + " " + str(xyz_insert[i][1]) + " " + 
           str(xyz_insert[i][2]))


   



sys.exit()

