#!/usr/bin/env python3
#
#    struc_insert: Inserts a molecule, cluster or surface
#      (given by POSCAR_insert) into a solvent, solid or 
#      similar structure (given by POSCAR_solvent) 
#      All atoms of the POSCAR_solvent located at the 
#      future position of POSCAR_insert will be removed.   
#    Part of VASP4CLINT
#     Julien Steffen, 2024 (julien.steffen@fau.de)
#

import sys 
import os 
import re
import numpy as np
from scipy.spatial import distance_matrix

print('''
This script inserts a given structure A (e.g., an intermetallic 
phase crystal) into another structure B (e.g., a SCALMS droplet).
For this, all atoms of B in the region where A is placed will be 
removed.
The unitcell of B should be much larger than A in order to contain it.
The dimensions of B will be taken as dimensions of the result.
Two POSCAR files are needed (both without selective dynamics!):
 - POSCAR_insert, containing the inserted structure A 
 - POSCAR_solvent, containing the solvent structure B 
The following keywords can/must be given:
 -center=x,y,z : Gives the center of mass position of the inserted 
    structure, within the solvent box. The values must be given in 
    Angstrom. Example: center=20.0,25.0,10.0
 -buffer=value : The distance between the solvent atoms and any atom
    of the inserted structure, below which the respective solvent 
    atom is deleted. Default: 2.8 (Angstrom). Example -buffer=3.5

The result will be written in POSCAR.
''')


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
#    the distance, below which atoms of the solvent are removed (in Angstrom)
#
dist_remove = 2.8

#
#    Always assume a global system scale of 1.0
#
sys_scale=1.0
#
#    Read in the command line arguments
#
read_center=False

for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-center":
         com_insert_read=actval.split(",")
         read_center=True
      if param == "-buffer":
         dist_remove=float(actval)
         bottom=True

if not read_center:
   print("Please give the desired insertion center of mass with the -center keyword!")
   sys.exit(1)

com_insert=np.zeros(3)
inc=0
for el in com_insert_read:
   com_insert[inc]=float(el)
   inc=inc+1

inc=inc-1
if inc < 2:
   print("Please give at least three values for the -center keyword!")
   sys.exit(1)
if inc > 2:
   print("Please give no more than three values for the -center keyword!")
   sys.exit(1)

print("Given settings:")
print(" - Inserted center of mass positions: x=",com_insert_read[0],", y=",com_insert_read[1],", z=",com_insert_read[2])
print(" - Buffer distance : ",str(dist_remove))
print(" ")
#
#    Read in the coordinates of structure A (inserted crystal or structure)
#

# Define coordinate transformation as functions in order to use them intermediately!
# F1: FRAC2CART
def trans_frac2cart(xyz,natoms,a_vec,b_vec,c_vec):
   xyz_cart = np.zeros((natoms,3))
   for i in range(natoms):
      xyz_cart[i][0]=(xyz[i][0]*a_vec[0]+xyz[i][1]*b_vec[0]+\
                     xyz[i][2]*c_vec[0])*sys_scale
      xyz_cart[i][1]=(xyz[i][0]*a_vec[1]+xyz[i][1]*b_vec[1]+\
                     xyz[i][2]*c_vec[1])*sys_scale
      xyz_cart[i][2]=(xyz[i][0]*a_vec[2]+xyz[i][1]*b_vec[2]+\
                     xyz[i][2]*c_vec[2])*sys_scale
            
   return xyz_cart
   
# F2: CART2FRAC
def trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec):
   xyz_frac = np.zeros((natoms,3))
   act_vec = np.zeros((3))
   vec_mat=np.zeros((3,3))
#
#   Build direct from cartesian coordinates by inverting the 
#   unit cell vector matrix and multiplying it to the coordinate vector
#

   vec_mat[0][0]=a_vec[0]
   vec_mat[1][0]=a_vec[1]
   vec_mat[2][0]=a_vec[2]
   vec_mat[0][1]=b_vec[0]
   vec_mat[1][1]=b_vec[1]
   vec_mat[2][1]=b_vec[2]
   vec_mat[0][2]=c_vec[0]
   vec_mat[1][2]=c_vec[1]
   vec_mat[2][2]=c_vec[2]

   mat_inv=np.linalg.inv(vec_mat)

   for i in range(natoms):
      xyz_frac[i][0]=(xyz[i][0]*mat_inv[0][0]+xyz[i][1]*mat_inv[0][1]+\
                     xyz[i][2]*mat_inv[0][2])*sys_scale
      xyz_frac[i][1]=(xyz[i][0]*mat_inv[1][0]+xyz[i][1]*mat_inv[1][1]+\
                     xyz[i][2]*mat_inv[1][2])*sys_scale
      xyz_frac[i][2]=(xyz[i][0]*mat_inv[2][0]+xyz[i][1]*mat_inv[2][1]+\
                     xyz[i][2]*mat_inv[2][2])*sys_scale


   return xyz_frac

print("Open file 'POSCAR_insert' ...")
insert_name="POSCAR_insert"

insert_in = open(insert_name,"r")

insert_cartesian=False
solvent_cartesian=False

with insert_in as infile:
   line = infile.readline()
   line = infile.readline()
   # read in the lattice vectors a,b and c

   a_insert=np.zeros(3)
   b_insert=np.zeros(3)
   c_insert=np.zeros(3)

   line = infile.readline().rstrip("\n")
   a_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   b_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   c_read = line.rstrip().split()[0:3]
   for i in range(3):
      a_insert[i]=float(a_read[i])
      b_insert[i]=float(b_read[i])
      c_insert[i]=float(c_read[i])

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
#    Detect if direct or cartesian coordinates are used
   if line == "Direct" or line == "direct":
      insert_cartesian=False
   else:
      insert_cartesian=True
 
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

#    If cartesian coordinates are used, transfer to direct
xyz_insert2 = np.zeros((natoms_insert,3)) 
if insert_cartesian:
   xyz_insert2=trans_cart2frac(xyz_insert,natoms_insert,a_insert,b_insert,c_insert)
else:
   xyz_insert2=xyz_insert
#    Move all atoms into unit cell 
for i in range(natoms_insert):
   for j in range(3):
      while xyz_insert2[i][j] < 0.0:
         xyz_insert2[i][j] = xyz_insert2[i][j] + 1.0
      while xyz_insert2[i][j] > 1.0:
         xyz_insert2[i][j] = xyz_insert2[i][j] - 1.0

#   Convert back to cartesian coordinates
xyz_insert=trans_frac2cart(xyz_insert2,natoms_insert,a_insert,b_insert,c_insert)

#
#    Read in the coordinates of structure B (solvent)
#

print("Open file 'POSCAR_solvent' ...")
solvent_name="POSCAR_solvent"

solvent_in = open(solvent_name,"r")

with solvent_in as infile:
   line = infile.readline()
   line = infile.readline()

   # read in the lattice vectors a,b and c

   a_solvent=np.zeros(3)
   b_solvent=np.zeros(3)
   c_solvent=np.zeros(3)

   line = infile.readline().rstrip("\n")
   unit_line1=line.rstrip("\n")
   a_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   unit_line2=line.rstrip("\n")
   b_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   unit_line3=line.rstrip("\n")
   c_read = line.rstrip().split()[0:3]
   for i in range(3):
      a_solvent[i]=float(a_read[i])
      b_solvent[i]=float(b_read[i])
      c_solvent[i]=float(c_read[i])

#   The elements line
   line = infile.readline()
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
      solvent_cartesian=False
   else:
      solvent_cartesian=True
   

   for i in range(natoms_solvent):
      line = infile.readline().rstrip("\n")
      xyz_read = line.rstrip().split()[0:3]
      for j in range(3):
         xyz_solvent[i][j]=float(xyz_read[j])

#    If cartesian coordinates are used, transfer to direct
xyz_solvent2 = np.zeros((natoms_solvent,3))
if solvent_cartesian:
   xyz_solvent2=trans_cart2frac(xyz_solvent,natoms_solvent,a_solvent,b_solvent,c_solvent)
else:
   xyz_solvent2=xyz_solvent
#    Move all atoms into unit cell 
for i in range(natoms_solvent):
   for j in range(3):
      while xyz_solvent2[i][j] < 0.0:
         xyz_solvent2[i][j] = xyz_solvent2[i][j] + 1.0
      while xyz_solvent2[i][j] > 1.0:
         xyz_solvent2[i][j] = xyz_solvent2[i][j] - 1.0

#   Convert back to cartesian coordinates
xyz_solvent=trans_frac2cart(xyz_solvent2,natoms_solvent,a_solvent,b_solvent,c_solvent)


#
#    Shift inserted structure to desired COM
#
com_shift=com_insert-com

for i in range(natoms_insert):
   for j in range(3):
      xyz_insert[i][j]=xyz_insert[i][j]+com_shift[j]


print("Compare solvent and insert structures, mark collisions ...")
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
print("Print final POSCAR file with combined structure ...")
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
   names_solv_print=[]
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
             names_solv_print.append(names_solvent[index])
             index2=index2+1
             elem_num_solv_new[i]=elem_num_solv_new[i]+1
             natoms_solv_new=natoms_solv_new+1
          index=index+1
#
#    Print the element names section
#    Combine both sections, each element shall only appear once
#

   elements_combined=[]
   for i in range (nelem_solvent):
      if elem_num_solv_new[i] > 0:
         elements_combined.append(elements_solvent[i])
   for i in range (nelem_insert):
      twice=False
      for el in elements_combined:
         if el == elements_insert[i]:
            twice=True
      if not twice:
         elements_combined.append(elements_insert[i])

   for el in elements_combined:
      print(el, end=" ")

   print(" ")

   


#
#    Print the element numbers section
#
   el_num_sum=[]
   for i in range(len(elements_combined)):
      el_num_sum.append(0)

   for i in range(len(elements_solvent)):
      for j in range(len(elements_combined)):
         if elements_combined[j] == elements_solvent[i]:
            el_num_sum[j]=el_num_sum[j]+elem_num_solv_new[i]

   for i in range(len(elements_insert)):
      for j in range(len(elements_combined)):
         if elements_combined[j] == elements_insert[i]:
            el_num_sum[j]=el_num_sum[j]+elem_num_insert[i]

   for ind in el_num_sum:
      print(ind, end=" ")
   print(" ")
 
   print("Cartesian")
   
   for ind in elements_combined:
  
      for i in range(natoms_solv_new):
         if names_solv_print[i] == ind:
            print(str(xyz_solv_print[i][0]) + " " + str(xyz_solv_print[i][1]) + " " + 
                 str(xyz_solv_print[i][2])) 

      for i in range(natoms_insert):
         if names_insert[i] == ind:
            print(str(xyz_insert[i][0]) + " " + str(xyz_insert[i][1]) + " " + 
                 str(xyz_insert[i][2]))

sys.stdout=original_stdout

print("Script struc_insert finished!")
print(" ")  



sys.exit()

