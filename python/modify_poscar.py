#!/usr/bin/env python3
#   
#    modify_poscar: Modify a given POSCAR file in several different
#      ways such as: shift its content, multiply it in different 
#      coordinate directions, convert it from frac 2 cart and the 
#      other way around, write it to a xyz file
#    Part of VASP4CLINT
#     Julien Steffen, 2023 (julien.steffen@fau.de)
#

import sys
import os
import re
import numpy as np
from numpy import linalg as LA

print('''
 This script takes a POSCAR file in direct or cartesian coordinates 
 and performs several different operations on it, depending 
 on the keyword and its specifiers given in the command line. 
 The list of possible options:
  -shift=a,b,c : Shift the unitcell contents along the given components 
     of a vector in direct coordinates along the coordinate axes.
     Example: shift=0.1,0.0,0.5 
  -multiply=a,b,c: Multiply the unit cell along a, b and c axes, where 
     the arguments are integers. Example: -multiply=2,2,1 (doubling the 
     cell and its content along x and y, keeping it along z)
  -frac2cart : If the POSCAR is in direct/fractional coordinates, convert
     it to cartesian coordinates
  -cart2frac : The the POSCAR is in cartesian coordinates, convert it to
     fractional coordinates 
  -writexyz : Write the coordinates of the POSCAR file to a xyz file  
 More than one job can be done at once, the ordering of operation is the 
 same as the ordering of keywords above 
''')

multiply_job=False
shift_job=False
frac2cart=False
cart2frac=False
writexyz=False

# Read in the command line arguments
for arg in sys.argv:
   if re.search("=",arg):
      param,actval=arg.split("=")
      if param == "-multiply":
         multiply_list=actval.split(",")
         multiply_job=True
      if param == "-shift":
         shift_list=actval.split(",")
         shift_job=True
   else:
      param=arg
      if param == "-frac2cart":
         frac2cart=True
      if param == "-cart2frac":
         cart2frac=True
      if param == "-writexyz":
         writexyz=True

if (not multiply_job) and (not shift_job) and (not frac2cart) and (not cart2frac) and (not writexyz):
   print("Please give a least one valid keyword!")
   sys.exit(1)

# Read in the POSCAR file
poscar_name="POSCAR"

poscar_in = open(poscar_name,"r")

# array for selective dynamics specifiers
coord_select= []  

cartesian=False
selective=False
with poscar_in as infile:
   line = infile.readline()
   line = line.rstrip("\n")
   poscar_comm = line # the comment line 
   line = infile.readline().rstrip("\n")
   sys_scale = float(line)  # the global lattice scaling factor
   # read in the lattice vectors a,b and c

   a_vec=np.zeros(3)
   b_vec=np.zeros(3)
   c_vec=np.zeros(3)

   line = infile.readline().rstrip("\n")
   a_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   b_read = line.rstrip().split()[0:3]
   line = infile.readline().rstrip("\n")
   c_read = line.rstrip().split()[0:3]
   for i in range(3):
      a_vec[i]=float(a_read[i])
      b_vec[i]=float(b_read[i])
      c_vec[i]=float(c_read[i])
   # read in the element ordering
   coord_vec=np.zeros(3)
   coord_vec[0]=a_vec[0]
   coord_vec[1]=b_vec[1]
   coord_vec[2]=c_vec[2]
   line = infile.readline().rstrip("\n")
   elements = line.rstrip().split()
   nelem = len(elements)
   # read in the number of elements 
   line = infile.readline().rstrip("\n")
   line_split = line.rstrip().split()

   elem_num=[]
   natoms=0
   names=[]
   for i in range(nelem):
      elem_num.append(int(line_split[i]))
      # total number of atoms
      natoms=natoms+elem_num[i]
      for j in range(elem_num[i]):
          names.append(elements[i])

   natoms=int(natoms)
   # read in the list of atoms 
   xyz = np.zeros((natoms,3))
   # check if selective dynamics was used
   line = infile.readline()
   line_split = line.rstrip().split()
   if (line_split[0] == "Selective" or line_split[0] == "selective"
        or line_split[0] == "Select" or line_split[0] == "select"):
      selective=True
      print(" The POSCAR file has selective dynamics.")
      print(" ")
   if selective:
      line = infile.readline().rstrip("\n")
   line_parts = line.split()
   if (line_parts[0] != "Direct" and line_parts[0] != "direct"  
                and line_parts[0] != " Direct" and line_parts[0] != " direct"):
      cartesian=True
      print(" The POSCAR has direct coordinates.")
   else:
      print(" The POSCAR has cartesian coordinates.")

   print(" ")
   for i in range(natoms):
      line = infile.readline().rstrip("\n")
      xyz_read = line.rstrip().split()[0:3]
      if selective:
         select_read = line.rstrip().split()[3:6]
         coord_select.append(select_read[0] + " " + select_read[1] + " " + select_read[2]) 
      for j in range(3):
         xyz[i][j]=float(xyz_read[j]) 

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


# A: SHIFT OF UNIT CELL ###############################
if shift_job:
   print(" Shift the unit cell along the vector: a=",shift_list[0],", b=",shift_list[1],
        " c=",shift_list[2],"  ...")
# iterating till the range
   shift_vec = []
   for i in range(3):
      ele = float(shift_list[i])
      shift_vec.append(ele) # adding the element



   if len(shift_vec) < 3:
      print("Please give three values for the coordinates!")
      exit(1)
   for i in range(3):
      if shift_vec[i] < -1.0:
         print("ERROR! All three shift coordinates must be greater than -1!")
         exit(1)
      if shift_vec[i] > 1.0:
         print("ERROR! All three shift coordinates must be smaller than 1!")
         exit(1)
   

   if selective:
      select_new=[]

   # Perform the actual shifting of the coordinates
   xyz_new=np.zeros((natoms,3))
   # If cartesian coordinates: translate the coordinates first to direct, then 
   #  perform the shift, and finally translate it back to cartesian!
   if cartesian:
      xyz = trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec)
      for i in range(natoms):
         for j in range(3):
            xyz_new[i][j] = xyz[i][j] + shift_vec[j]
            if xyz[i][j] < 0.0:
               xyz_new[i][j] = xyz_new[i][j] + 1.0
            if xyz[i][j] > 1.0:
               xyz_new[i][j] = xyz_new[i][j] - 1.0
         if selective:
            select_new.append(coord_select[i])
      xyz_new = trans_frac2cart(xyz_new,natoms,a_vec,b_vec,c_vec)

   else:
      for i in range(natoms):
         for j in range(3):
            xyz_new[i][j] = xyz[i][j] + shift_vec[j]
            if xyz[i][j] < 0.0:
               xyz_new[i][j] = xyz_new[i][j] + 1.0
            if xyz[i][j] > 1.0:
               xyz_new[i][j] = xyz_new[i][j] - 1.0 

         if selective:
            select_new.append(coord_select[i])  
   print(" done!\n")

# B: MULTIPLY THE UNITCELL #################################
if multiply_job:
   print(" Multiply the unit cell: ",multiply_list[0],"x along a,",multiply_list[1],
                 "x along b,",multiply_list[2],"x along c ...") 
   mult_vec = np.zeros(3)
   for i in range(3):
      ele = int(multiply_list[i])
      mult_vec[i] = ele # adding the element
#   If the previous job (shift of unit cell) was already done, overwrite the xyz array   
   if shift_job:
      xyz=xyz_new

#   First, determine the size and the number of atoms in the new unit cell      
  
   factor=int(mult_vec[0])*int(mult_vec[1])*int(mult_vec[2])
   natoms_old=natoms
   natoms=natoms*factor

   for i in range(nelem):
      elem_num[i]=elem_num[i]*factor


#  Scale the atom coordinates
   xyz_new=np.zeros((natoms,3))
#  The individual atoms: their elements   
   names_new=[]
   if selective:
      select_new=[]
   if cartesian:
      pos_new=0
      pos_old=0
      for i in range(natoms_old):
         for j in range(int(mult_vec[0])):
            for k in range(int(mult_vec[1])):
                for l in range(int(mult_vec[2])):
                   xyz_new[pos_new,0]=(xyz[pos_old,0]+float(j)*a_vec[0]+float(k)*b_vec[0]+
                                       float(l)*c_vec[0])
                   xyz_new[pos_new,1]=(xyz[pos_old,1]+float(j)*a_vec[1]+float(k)*b_vec[1]+
                                       float(l)*c_vec[2])
                   xyz_new[pos_new,2]=(xyz[pos_old,2]+float(j)*a_vec[2]+float(k)*b_vec[2]+
                                       float(l)*c_vec[2])
                   pos_new=pos_new+1
                   if selective:
                      select_new.append(coord_select[pos_old])
                   names_new.append(names[pos_old])   
         pos_old=pos_old+1


   else:
      pos_new=0
      pos_old=0
      for i in range(natoms_old):
         for j in range(int(mult_vec[0])):
            for k in range(int(mult_vec[1])):
                for l in range(int(mult_vec[2])):
                   xyz_new[pos_new,0]=(xyz[pos_old,0]+1.0*float(j))/float(mult_vec[0])
                   xyz_new[pos_new,1]=(xyz[pos_old,1]+1.0*float(k))/float(mult_vec[1])
                   xyz_new[pos_new,2]=(xyz[pos_old,2]+1.0*float(l))/float(mult_vec[2])
                   pos_new=pos_new+1
                   if selective:
                      select_new.append(coord_select[pos_old])
                   names_new.append(names[pos_old])
         pos_old=pos_old+1

#   Scale the unit cell vectors
   a_vec=np.multiply(a_vec,mult_vec[0])
   b_vec=np.multiply(b_vec,mult_vec[1])
   c_vec=np.multiply(c_vec,mult_vec[2])

   names=names_new
   coord_select=select_new
   print(" done!\n")
# C: TRANSLATE FROM DIRECT TO CARTESIAN ###################################
if frac2cart:
   print(" Translate struture from fractional/direct to cartesian coordinates...")   
   if cartesian:
      print(" The POSCAR file already has cartesian coordinates!")
      exit(1)
#  If one of the previous jobs were done, overwrite the xyz array
   if shift_job or multiply_job:
      xyz=xyz_new 

   xyz_new=trans_frac2cart(xyz,natoms,a_vec,b_vec,c_vec)
   cartesian=True
   select_new=[]
   if selective:
      for i in range(natoms):
         select_new.append(coord_select[i])


   print(" done!\n")
  
# D: TRANSLATE FROM CARTESIAN TO DIRECT ###################################
if cart2frac:
   print(" Translate struture from cartesian to fractional/direct coordinates...")
   if not cartesian:
      print(" The POSCAR file already has direct coordinates!")
      exit(1)    
#  If one of the previous jobs were done, overwrite the xyz array
   if shift_job or multiply_job or frac2cart:
      xyz=xyz_new


   xyz_new=trans_cart2frac(xyz,natoms,a_vec,b_vec,c_vec)

   cartesian=False
   select_new=[]
   if selective:
      for i in range(natoms):
         select_new.append(coord_select[i])


   print(" done!\n")

# E: WRITE XYZ FILE #######################################################
if writexyz:
   print(" Write structure to xyz file poscar_mod.xyz...")
#  If no other job has been done before, copy xyz directly to xyz_new
   if not shift_job and not multiply_job and not frac2cart and not cart2frac:
      xyz_new=xyz 
#  If the structure is in direct coordinates, translate it first to cartesian!
   if (not cartesian):
      xyz_new=trans_frac2cart(xyz_new,natoms,a_vec,b_vec,c_vec) 

   original_stdout=sys.stdout
   with open("poscar_mod.xyz","w") as f:
      sys.stdout = f
      print(natoms)
      print(" System converted to xyz by modify_poscar.py")
      for i in range(natoms):
         print(names[i],"  ",str(xyz_new[i][0]),"  ",str(xyz_new[i][1]),"  ",
                    str(xyz_new[i][2]))
#  Translate the structure back to direct coordinates, for final POSCAR printout          
   if (not cartesian):
      xyz_new=trans_cart2frac(xyz_new,natoms,a_vec,b_vec,c_vec)

   sys.stdout=original_stdout
   print(" done!\n")
    




 
# Finally, write out the modified POSCAR file POSCAR_mod

print(" New file POSCAR_mod will be written...")

original_stdout=sys.stdout
with open("POSCAR_mod","w") as f:
   sys.stdout = f

   # rescale the z-component of the unit cell vector such that the desired 
   # total length in z-direction of the system is retained

   print("Modified POSCAR written by modify_poscar.py")
   print(sys_scale)
   print(str(a_vec[0]) + "  " + str(a_vec[1]) + "  " + str(a_vec[2]))
   print(str(b_vec[0]) + "  " + str(b_vec[1]) + "  " + str(b_vec[2]))
   print(str(c_vec[0]) + "  " + str(c_vec[1]) + "  " + str(c_vec[2]))
   for i in range(nelem):
      f.write(elements[i] + " ")
   print(" ")
# Combine element sections of both molecules  
   elem_num_all=[]
   elements_all =[]
   for i in range(nelem):
      elem_num_all.append(elem_num[i])
   for elem in elem_num_all:
      f.write(str(elem) + " ")
   
   print(" ")
# Write combined element numbers section
   if selective == 1:
      print("Selective Dynamics")
   if cartesian:
      print("Cartesian")
   else:
      print("Direct")

   if selective:
      for i in range(natoms):
         print(str(xyz_new[i][0]) + " " + str(xyz_new[i][1]) + " " +
              str(xyz_new[i][2]) + " " + select_new[i])
   else:
      for i in range(natoms):
         print(str(xyz_new[i][0]) + " " + str(xyz_new[i][1]) + " " +
              str(xyz_new[i][2])) 

sys.stdout=original_stdout
print(" done!")

print ('''
 Execution of modify_poscar.py successfully finished!
 File POSCAR_mod has been written!
   ''')
