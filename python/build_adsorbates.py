#!/usr/bin/env python3
#   
#    build_adsorbates: Place adsorbate molecules or atoms 
#      on surfaces given by a POSCAR file. The adsorbates can be 
#      chosen from a predefined library or given locally 
#    Part of VASP4CLINT
#     Julien Steffen, 2023 (julien.steffen@fau.de)
#

import sys
import os
import numpy as np
import random as ran
import io
from numpy import linalg as LA
from scipy.spatial import distance_matrix

# This script manages the complete builtup of VASP
# input for molecules adsorbed on a surface systems.
# The POSCAR file of a metal surface slab and the POSCAR
# or XYZ files of the ionic liquids on it (anion and cation,
# respectively) need to be provided, as well as the location
# of the directory containing the potcar directory tree 
# from VASP.

print('''
 This script manages the complete buildup of VASP
 input for systems with molecules adsorbed on surfaces.
 The POSCAR file of a (metal) surface slab (POSCAR_surf) needs
 to be provided. POSCAR files of adsorbed molecules should
 be located in a folder with suitable names (look into the
 script header). A directory for the POTCAR files for the
 elements must be given there as well.
 Add -select_dist, if only metal atoms near to the molecule shall
   be optimized (change opt_dist parameter in script)
 Add -rot_surf, if the metal surface shall be rotated by 90 degrees
 Add -random, if the rotations shall be chosen by chance and also the
   positions shall be varied by +/- 0.5 Angstrom in x- and y-direction
 Add -cutout, if atoms of the surface near enough to any adsorbate
   shall be removed to avoid collisions (e.g. for placement in surface)
   (works only for adsorbates in principial unit cell!)
 Further keywords for submission of individual values:
 -z_length=[value]: Desired vacuum along z-axis in Angstrom (default: 23)
 -dist_cut=[value]: Desired distance to adsorbate atoms below which atoms
   of the surface are removed for being too near (for -cutout active)
 The list of molecules will be read in from file 'adsorb_list.dat'
 The format shall be: [Species No or shortcut] [x-coord] [y-coord] 
 [dist mol-surf(z)] [prec.] [nut.] [rot.] (angles given in degrees)"
''')


# GENERAL PARAMETERS:

# If select_dist is active: the maximum distance of atoms free for optimization
# to the adatoms
opt_dist = 7.0

# The absolute path to the ionic liquid molecules
mol_path="/home/jsteffen/work/build_adsorbates_input/"

# The path where the POTCAR files for the elements are located
potcar_path = "/scratch/potcar/PAW_PBE.52/"

# Total number of adsorbates available (besides individual molecule)
adsorb_num = 8

# lists of available ionic liquids
mol_names = [] 
mol_names.append("OTf (Triflate anion)")
mol_names.append("EMIM (Ethylmethylimidazolium cation)")
mol_names.append("NTF2 (Bistriflimide anion)")
mol_names.append("C4C1Pyr (Butylmethylimidazolium cation)")
mol_names.append("CO (Cabonmomoxide)")
mol_names.append("H2 (Hydrogen molecule)")
mol_names.append("H1 (Hydrogen atom)")
mol_names.append("Ga1 (Gallium atom)")
mol_names.append("POSCAR_ads from local folder")

mol_short = []
mol_short.append("otf")
mol_short.append("emim")
mol_short.append("ntf2")
mol_short.append("c4c1pyr")
mol_short.append("co")
mol_short.append("h2")
mol_short.append("h1")
mol_short.append("ga1")
mol_short.append("local")

# POSCAR files with adsorbate molecules
adsorbate_dict = {'otf' : '''Shifted unit cell system
1.0
10.0  0.0  0.0
0.0  10.0  0.0
0.0  0.0  10.0
S C O F
1 1 3 3
Direct
0.527015183558875 0.4607082659323616 0.31784703203214737
0.524609166818558 0.6439381998229663 0.2732994187024371
0.5406443096943152 0.3972853431046388 0.18718855081328112
0.6437205878251322 0.4516607494438104 0.4048892858087037
0.3978655913836207 0.44318043979277455 0.38333654006595697
0.6450714409634255 0.6868610719056076 0.2242203689640525
0.4959073480187655 0.72302616149723 0.3813361505269446
0.4308843717373069 0.6727537685006032 0.1776506530864752''',\
        'emim': '''Shifted unit cell system
1.0
20.0  0.0  0.0
0.0  20.0  0.0
0.0  0.0  20.0
H C N
13 6 2
Direct
0.25541869336606704 0.5071182614447782 0.5756952035706403
0.28596751209676263 0.5270220536744341 0.6573261830056848
0.33079527949382026 0.6256346235816651 0.6026677752452889
0.3297185936319754 0.588175342093485 0.5228418279886918
0.45710385889519356 0.46201694495909074 0.5878371929528976
0.4520497557122535 0.6101806734492015 0.6633489725523879
0.4555047597927526 0.6533319168986883 0.5859704788617388
0.5384961056810762 0.5258456784919204 0.6073571163745477
0.5676026796415556 0.6036632728251126 0.633814559749385
0.5510171281815437 0.5896408277139608 0.5479229007726122
0.28997317650091836 0.37839611818736696 0.5915267566226933
0.3450925509888232 0.3777740201396833 0.6621110558819041
0.3766404689646037 0.3629257749053748 0.5802783456197623
0.300472484374958 0.5192524127213536 0.6045009652688513
0.33986004905698064 0.5790846234248038 0.5761082140683664
0.41121614652643734 0.4917244753225436 0.5906673228412961
0.46261802514548844 0.6042772264366134 0.6096862894804579
0.5334649418178627 0.5795973225248524 0.5986125747517878
0.3400218229444465 0.39200631559932764 0.6090216376816193
0.351440256283271 0.4640712527198715 0.5988616396927227
0.41057471090321784 0.5577968628858657 0.5831964870166698''',\
        'ntf2': '''Converted_POSCAR
1.0
+10.0000000000  +0.0000000000  +0.0000000000 
 +0.0000000000 +10.0000000000  +0.0000000000
 +0.0000000000  +0.0000000000 +10.0000000000
N O F C S
 1 4 6 2 2
Cartesian
 -0.6916073607  +1.2731033737  -0.2234933979
 -2.0257149253  +2.1854328958  +1.6404658470
 -3.1351199623  +1.0945274766  -0.3214276384
 +1.4388033365  +0.9944187336  +1.0709925481
 +0.6560407781  -0.8737283678  -0.3943191261
 -1.0014130834  -0.6350227309  +2.3293310497
 -2.1282167516  -1.5198313496  +0.6568264536
 -3.1999209207  -0.5224358809  +2.2981484056
 +2.8540310209  +1.0015672864  -1.5861494394
 +1.4487476677  +2.6938871102  -1.5318182856
 +0.9788395518  +0.9139410293  -2.7362563685
 -2.0867119440  -0.5114546904  +1.5420266476
 +1.5545319418  +1.3570404913  -1.6049956141
 +0.7379110758  +0.5381487139  -0.1078039484
 -2.0537204246  +1.1546759088  +0.6297128669''',\
         'c4c1pyr': '''Converted_POSCAR
   1.00000000000000
    20.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   20.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   20.0000000000000000
   N    C    H
     1     9    20
Direct
  0.3296653713507665  0.5958078884348149  0.4748624124652654
  0.3658327565996023  0.6612240744585941  0.4929627898719906
  0.3157186149127540  0.7168475063746764  0.4793513488909960
  0.2490893981190105  0.6843196792136375  0.4986305010505493
  0.2551400426245072  0.6141542945817203  0.4696128952485760
  0.3399260557863408  0.5451600093279039  0.5293927423086443
  0.3558880897667929  0.5654775470344240  0.4100165904176441
  0.3587785304356891  0.6119008988898738  0.3497677212024224
  0.3867158915193291  0.5743557848544774  0.2888931346689793
  0.3915318923222632  0.6193657271336727  0.2273092833648514
  0.4130923003091197  0.6636176906195587  0.4648617057075083
  0.3778037194685530  0.6583256251957748  0.5465006002545176
  0.3157661917494520  0.7313461997986264  0.4263043011635388
  0.3275878104815886  0.7614462688482417  0.5092515394860844
  0.2051575585917330  0.7102868292261439  0.4780529052320374
  0.2437417699961886  0.6823823323970678  0.5534253909992370
  0.2415936117886018  0.6130502740461174  0.4164189824720221
  0.2260788782991969  0.5754311294241911  0.4957292420513325
  0.3940235024528672  0.5372008568149905  0.5360372724204093
  0.3151422082594764  0.4982615723704837  0.5147969750155207
  0.3177760861645558  0.5644423593294853  0.5758162622284558
  0.3233009907745145  0.5221192050277889  0.3999162652385450
  0.4062631742459456  0.5466746801732412  0.4219644432970935
  0.3911249919220048  0.6552323948848051  0.3601956221832863
  0.3088243549780276  0.6315956840490702  0.3374090338416971
  0.3546115759575239  0.5308796828731123  0.2779546650992656
  0.4365430010641940  0.5542415825906668  0.3012955646054259
  0.3421524526385581  0.6391803344153900  0.2132888452616809
  0.4114274611145116  0.5918770037247538  0.1841089901207100
  0.4247587163063396  0.6622258838867113  0.2366384738317095''',\
        'co' : ''' fcc (111) surface
   1.00000000000000     
    20.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   20.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   20.0000000000000000
   C    O
     1     1
Direct
  0.0688774879481309  0.0791554878970402  0.5262517800214072
  0.0687629486018699  0.0795448747879633  0.5834259439785937''',\
          'h2' : '''Converted_POSCAR
   1.00000000000000
    20.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   20.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   20.0000000000000000
   H
     2
Direct
  0.5039870221033285  0.5004552151826146  0.5000447848173878
  0.5415129778966685  0.5000447848173878  0.5004552151826146''',\
          'h1' : '''Converted_POSCAR
   1.00000000000000
    20.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   20.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   20.0000000000000000
   H
     1
Direct
  0.5415129778966685  0.5000447848173878  0.5004552151826146''',\
          'ga1' : '''Converted_POSCAR
   1.00000000000000
    20.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   20.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   20.0000000000000000
   Ga
     1
Direct
  0.5415129778966685  0.5000447848173878  0.5004552151826146'''}





# Molecular weights of ionic liquid molecules for bulk buildup
mol_weights = np.zeros(2)
mol_weights[0] = 149.0632
mol_weights[1] = 113.184


# The elongation of the system in z-direction
z_length=23.0

# The distance between the molecule and the surface
mol_surf_dist=2.3

# distance for cutout of surface atoms, if the option "cutout" is activated
dist_remove=1.8

# Check if command line arguments were passed

# If needed, activate the flag that only z-coordinates shall be 
# optimized
zonly = False
bulk = False
random = False
rot_surf = False
cutout = False
select_dist = False
if len(sys.argv) > 1:
   for arg in sys.argv: 
      if arg == "-zonly":
         zonly = True 
         print("Only z coordinates will be optimized.")
# If the metal surface shall be rotated by 90 degrees for detection 
# of herringbone patterns
      if arg == "-rot_surf":
         rot_surf = True 
         print("The underlying surface will be rotated by 90 degrees.")
# If the positions of the molecules on the surface shall be determined
# (partially) by chance. The rotations will be chosen as random numbers,
# the given translations in the file will be varied by up to +/- 0.5 A,
      if arg == "-random":
         random = True
         print("The positions of the molecules on the surface shall be determined") 
         print("by chance. Rorations random, translation given values +/. 0.5 A.")
# If atoms of the support too near to any adsorbate shall be cut out
      if arg == "-cutout":
         cutout = True
         print("Surface atoms too near to adsorbate will be cut out.")
# Switch for selective dynamics based on distance to additive atoms 
      if arg == "-select_dist":
          select_dist = True
          print("Only atoms near enough to adsorbates will be free for optimization.")
# Read in the desired length of the system along z-axis
      if arg[0:10] == "-z_length=":
          z_length = float(arg[10:])
          print("Length of vacuum in system along z-axis changed to " + str(z_length) + "A.")
# Read in the distance to adsorbates below which surface atoms will be removed
# (only if -cutout option is activated!)
      if arg[0:10] == "-dist_cut=":
          dist_remove = float(arg[10:])
          print("Distance for surface atom removal changed to " + str(dist_remove) + " A.")


# Predefine the writing POSCAR flags for selective dynamics
if zonly:
   select_true = " F F T "
else:
   select_true = " T T T "
select_false = " F F F "


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
#   Calculate Euler rotation matrix for given angles (precision, nutation, rotation)
#
def EulerMatrix(pre,nut,rot):
   Emat=np.zeros((3,3))
   Emat[0][0] = np.cos(rot)*np.cos(pre)-np.sin(rot)*np.cos(nut)*np.sin(pre)
   Emat[0][1] = np.cos(rot)*np.sin(pre)+np.sin(rot)*np.cos(nut)*np.cos(pre)
   Emat[0][2] = np.sin(rot)*np.sin(nut)
   Emat[1][0] = -np.sin(rot)*np.cos(pre)-np.cos(rot)*np.cos(nut)*np.sin(pre)
   Emat[1][1] = -np.sin(rot)*np.sin(pre)+np.cos(rot)*np.cos(nut)*np.cos(pre)
   Emat[1][2] = np.cos(rot)*np.sin(nut)
   Emat[2][0] = np.sin(nut)*np.sin(pre)
   Emat[2][1] = -np.sin(nut)*np.cos(pre)
   Emat[2][2] = np.cos(nut)
   return Emat

#
#   Calculate the center of Mass of the surface
#
def COM(xyz,natoms,names,elements_dict,move):
   com=np.zeros(3)
   mass=0.0
   for i in range(natoms):
      for j in range(3):
         com[j]=com[j]+elements_dict[names[i]]*xyz[i][j]
      mass=mass+elements_dict[names[i]]

   for i in range(3):
      com[i]=com[i]/mass
   if (move):
      for i in range(natoms):
         for j in range(3):
            xyz[i][j]=xyz[i][j]-com[j]
   return com      

#
#   Rotate a molecule with respect to the given rotation matrix
#
def RotMolecule(xyz,natoms,matrix):
   xyz_rot = np.zeros((natoms,3))
   for i in range(natoms):
      coord = np.zeros(3)
      coord_new = np.zeros(3)
      for j in range(3):
         coord[j]=xyz[i][j]
         coord_new=np.dot(matrix,coord)
      for j in range(3):
         xyz_rot[i][j]=coord_new[j]
   return xyz_rot
#
#   Print a molecule xyz to the file with given name
#   Parameters:
#   - xyz: the xyz coordinates of the molecule
#   - names: the element symbols of the atoms
#   - natoms: number of atoms in the molecule 
#   - comment: Desired entries in the xyz comment line
#   - filename: name of output file (xyz)
#
def PrintMolecule(xyz,names,natoms,comment,filename):
   original_stdout=sys.stdout
   with open(filename,"w") as f:
      sys.stdout = f
      print(natoms)
      print(comment)
      for i in range(natoms):
         print(names[i] + "  " +  str(xyz[i][0]) + "  " + str(xyz[i][1])
                  + "  " + str(xyz[i][2]))

      sys.stdout=original_stdout
#############################################################################
#    Class molecule: for ionic liquid molecules on the surface 
#
class IL_Molecule:
#
#    The general constructor
#
   def __init__(self, IL_spec, IL_xpos0, IL_ypos0, IL_surf_dist,IL_prec, IL_nut, IL_rot):
      self.spec = IL_spec
      self.xpos0 = IL_xpos0
      self.ypos0 = IL_ypos0
      self.surf_dist = IL_surf_dist
      self.prec = IL_prec
      self.nut = IL_nut
      self.rot = IL_rot
#
#    Read in the structure of the molecule
#
   def ReadMolecule(self, mol_name):

#
#    Either read from the local file POSCAR ads or from the dictionary of the 
#    python script itself!
#
    
      if mol_name == "POSCAR_ads":
         mol_in = open(mol_name,"r") 
      elif mol_name == "POSCAR_ads2":
         mol_in = open(mol_name,"r") 
      else:   
         mol_in = io.StringIO(adsorbate_dict[mol_name])

      with mol_in as infile:
         line = infile.readline().rstrip("\n")

         line = infile.readline().rstrip("\n")
         mol_scale = float(line)  # the global lattice scaling factor
   # read in the lattice sizes (cubic primitive assumed)

         a_vecm=np.zeros(3)
         b_vecm=np.zeros(3)
         c_vecm=np.zeros(3)

         line = infile.readline().rstrip("\n")
         a_read = line.rstrip().split()[0:3]
         line = infile.readline().rstrip("\n")
         b_read = line.rstrip().split()[0:3]
         line = infile.readline().rstrip("\n")
         c_read = line.rstrip().split()[0:3]
         for i in range(3):
            a_vecm[i]=float(a_read[i])
            b_vecm[i]=float(b_read[i])
            c_vecm[i]=float(c_read[i])
         xlen=a_vecm[0]
         ylen=b_vecm[1]
         zlen=c_vecm[2]
   # read in the element ordering
         line = infile.readline().rstrip("\n")
         self.elements = line.rstrip().split()
         self.nelem = len(self.elements)
   # read in the number of elements 
         line = infile.readline().rstrip("\n")
         line_split = line.rstrip().split()

         self.elem_num=[]
         self.names=[]  # array with element symbols of the molecule (testing)
         self.natoms=0
         for i in range(self.nelem):
            self.elem_num.append(int(line_split[i]))
            for j in range(self.elem_num[i]):
                self.names.append(self.elements[i])
      # total number of atoms in the surface
            self.natoms=self.natoms+self.elem_num[i]

   # read in the list of atoms in the first molecule
         self.xyz = np.zeros((self.natoms,3))
   # one line less, since no selective dynamics is needed!
   # if a index out of range error occures here, change the POSCAR file!
         line = infile.readline()
         line_split = line.rstrip().split()
   # Check if direct or cartesian coordinates are used!
         coord_direct = True
         if line_split[0] == "Cartesian":
            coord_direct = False
         elif line_split[0] == "Direct":
            coord_direct = True
         for i in range(self.natoms):
            line = infile.readline().rstrip("\n")
            xyz_read = line.rstrip().split()[0:3]
            for j in range(3):
               self.xyz[i][j]=float(xyz_read[j])
            if coord_direct:
               self.xyz[i][0]=self.xyz[i][0]*xlen
               self.xyz[i][1]=self.xyz[i][1]*ylen
               self.xyz[i][2]=self.xyz[i][2]*zlen


   # Remove periodicity from atoms, move molecule near to origin
      if coord_direct:
         
         for i in range(self.natoms):
            if self.xyz[i][0] > xlen*0.78:
               self.xyz[i][0]=(self.xyz[i][0]-xlen)
            if self.xyz[i][1] > ylen*0.78:
               self.xyz[i][1]=(self.xyz[i][1]-ylen)
            if self.xyz[i][2] > zlen*0.78:
               self.xyz[i][2]=(self.xyz[i][2]-zlen)

#
#   Calculate the center of mass (COM) of a molecule
#   Optional flag: move the molecule to the COM if desired
#
   def COM_Molecule(self):
      self.com=np.zeros(3)
      self.mass=0.0
      for i in range(self.natoms):
         for j in range(3):
            self.com[j]=self.com[j]+elements_dict[self.names[i]]*self.xyz[i][j]
         self.mass=self.mass+elements_dict[self.names[i]]

      for i in range(3):
         self.com[i]=self.com[i]/self.mass
      if (True):
         for i in range(self.natoms):
            for j in range(3):
               self.xyz[i][j]=self.xyz[i][j]-self.com[j]

#
#   Calculate the moments of inertia tensor of the molecule 
#
   def MIT_Molecule(self):

      self.MIT=np.zeros((3,3))
      for i in range(self.natoms):
   # the diagonal elements
         self.MIT[0][0]=self.MIT[0][0]+elements_dict[self.names[i]]*\
                  (self.xyz[i][1]**2+self.xyz[i][2]**2)
         self.MIT[1][1]=self.MIT[1][1]+elements_dict[self.names[i]]*\
                  (self.xyz[i][0]**2+self.xyz[i][2]**2)
         self.MIT[2][2]=self.MIT[2][2]+elements_dict[self.names[i]]*\
                  (self.xyz[i][0]**2+self.xyz[i][1]**2)
   # the off-diagonal elements
         self.MIT[0][1]=self.MIT[0][1]-elements_dict[self.names[i]]*\
                  (self.xyz[i][0]*self.xyz[i][1])
         self.MIT[0][2]=self.MIT[0][2]-elements_dict[self.names[i]]*\
                  (self.xyz[i][0]*self.xyz[i][2])
         self.MIT[1][2]=self.MIT[1][2]-elements_dict[self.names[i]]*\
                  (self.xyz[i][1]*self.xyz[i][2])
# the mirrored off-diagonal elements
         self.MIT[1][0]=self.MIT[0][1]
         self.MIT[2][0]=self.MIT[0][2]
         self.MIT[2][1]=self.MIT[1][2]

#
#   Rotate a molecule with respect to the given rotation matrix
#
   def Rot_Molecule(self,matrix):
      xyz_rot = np.zeros((self.natoms,3))
      for i in range(self.natoms):
         coord = np.zeros(3)
         coord_new = np.zeros(3)
         for j in range(3):
            coord[j]=self.xyz[i][j]
            coord_new=np.dot(matrix,coord)
         for j in range(3):
            xyz_rot[i][j]=coord_new[j]
      self.xyz = xyz_rot
#  End of Class molecule
######################################################################


# First, read in the POSCAR file for the metal surface
surface_name="POSCAR_surf"

surface_in = open(surface_name,"r")

with surface_in as infile:
   line = infile.readline()
   line = line.rstrip("\n")
   surf_comm = line # the comment line 
   line = infile.readline().rstrip("\n")
   surf_scale = float(line)  # the global lattice scaling factor
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

   # calculate normalized vectors for bulk builtup
   # take absolute values for definition of border planes
   a_len = np.linalg.norm(a_vec)
   b_len = np.linalg.norm(b_vec)
   c_len = np.linalg.norm(c_vec)
   a_norm = np.zeros(3)
   b_norm = np.zeros(3)
   c_norm = np.zeros(3)
   a_norm[0] = abs(a_vec[1])/a_len
   a_norm[1] = abs(a_vec[0])/a_len
   a_norm[2] = abs(a_vec[2])/a_len
   b_norm[0] = abs(b_vec[1])/b_len
   b_norm[1] = abs(b_vec[0])/b_len
   b_norm[2] = abs(b_vec[2])/b_len
   c_norm[0] = abs(c_vec[0])/c_len
   c_norm[1] = abs(c_vec[1])/c_len
   c_norm[2] = abs(c_vec[2])/c_len

   

   # read in the element ordering
   line = infile.readline().rstrip("\n")
   elements_surf = line.rstrip().split()
   nelem = len(elements_surf)
   # read in the number of elements 
   line = infile.readline().rstrip("\n")
   line_split = line.rstrip().split()

   elem_num_surf=[]
   natoms_surf=0
   names_surf=[]
   for i in range(nelem):
      elem_num_surf.append(int(line_split[i]))
      # total number of atoms in the surface
      natoms_surf=natoms_surf+elem_num_surf[i]
      for j in range(elem_num_surf[i]):
          names_surf.append(elements_surf[i])

   natoms_surf=int(natoms_surf)
   # read in the list of atoms in the surface 
   xyz_surf = np.zeros((natoms_surf,3))
   line = infile.readline()
   line=line.rstrip()
   line_split=line.split()
   if line_split[0] == "Selective" or line_split[0] == "selective":
      line = infile.readline()

   line=line.rstrip()
   if line_split[0] == "Direct" or line_split[0] == "direct":
      print(" Please give POSCAR_surf in cartesian coordinates!")
      sys.exit()    
   for i in range(natoms_surf):
      line = infile.readline().rstrip("\n")
      xyz_read = line.rstrip().split()[0:3]
      for j in range(3):
         xyz_surf[i][j]=float(xyz_read[j])

# Now read in the structure of ionic liquid molecules 
# The user can select from a predefined library of molecules for that purpose
   
print(" List of available adsorbate molecules (with No. and shortcut):")
for i in range (adsorb_num):
    print(" -(" + str(i+1) +" / " + mol_short[i] + ")  --  "+ mol_names[i])

print(" ")
print(" -(99)- Individual molecule (POSCAR_ads in current folder)")
print(" -(100)- Individual molecule (POSCAR_ads2 in current folder)")
print(" ")

nspecs=0
#while True:
#   try: 
#      nspecs = int(input(" Number of ionic liquid molecules?  "))
#      if nspecs <= 0 or nspecs > 8:
#          raise ValueError
#      break 
#   except NameError:
#      print(" ERROR! No valid number given! Please try again!")
#   except ValueError:
#      print(" ERROR! Only numbers between 1 and 8 are allowed!")

# Read in the species and positions/rotations from file!

# Second, read in the data file with the bromine positions
list_name="adsorb_list.dat"

list_in = open(list_name,"r")

with list_in as infile:
   IL_lines = list_in.readlines()

nspecs = len(IL_lines)
#if len(IL_lines) != nspecs:
#   print("Error! The file IL_list.dat must contain as many lines as the")
#   print("  number of desired IL molecules.")

#
#    Now read in all surface molecules as class instances!
#    First, only their general specifiers
#

all_mols = []
for line in (IL_lines):
   pos_read = line.rstrip().split()[0:7]
#
#    Both read in variants: direct index or over name shortcut
#
   try:
      IL_spec = (int(pos_read[0]))
   except ValueError:
      IL_spec = mol_short.index(str(pos_read[0]))+1 

   if IL_spec <= 0 or IL_spec > adsorb_num and IL_spec != 99 and IL_spec != 100:
      print("Error! Only indices between 1 and " + str(adsorb_num) + " (or 99 and 100) allowed!")    

   IL_xpos0 = (float(pos_read[1]))
   IL_ypos0 = (float(pos_read[2]))
   IL_surf_dist = (float(pos_read[3]))
   IL_prec = (float(pos_read[4])/180.0*np.pi)
   IL_nut = (float(pos_read[5])/180.0*np.pi)
   IL_rot = (float(pos_read[6])/180.0*np.pi)

#
#    If random mode is activated, determine angles by chance and vary x and y 
#    positions partially by chance (+/- 0.5 A)
#
   if random:
      IL_prec = ran.uniform(0.0, 360.0)/180.0*np.pi 
      IL_nut = ran.uniform(0.0, 180.0)/180.0*np.pi  # range from 0 to 180!
      IL_rot = ran.uniform(0.0, 360.0)/180.0*np.pi
      IL_xpos0 = IL_xpos0 + ran.uniform(-0.5, 0.5)
      IL_ypos0 = IL_ypos0 + ran.uniform(-0.5, 0.5)

   mol_act = IL_Molecule(IL_spec, IL_xpos0, IL_ypos0, IL_surf_dist, IL_prec, IL_nut, IL_rot)

   all_mols.append(mol_act)

#  if the ionic liquid bulk shall be placed above the surface, being 
if 1==2:
   print("bullshit!")
else :

#
#    Loop through all molecules in the list and determine their structures!
#

   for m in range(nspecs): 
# Read in the first molecule from its POSCAR in the build_scill_input folder
   
# First, read in the POSCAR file for the metal surface
      if all_mols[m].spec-1 == 98:
         mol_name="POSCAR_ads"
      elif all_mols[m].spec-1 == 99:
         mol_name="POSCAR_ads2" 
      else:   
         mol_name=mol_short[all_mols[m].spec-1]

# Call subroutine for read in of the molecule
      all_mols[m].ReadMolecule(mol_name)


# Now rotate molecule to inertia axis frame, in order to be capable 
# of depositing it on the surface
# Determine center of masss
      all_mols[m].COM_Molecule()

# Calculate the moments of inertia tensor of the molecule 

      all_mols[m].MIT_Molecule()
# solve the eigenvalue problem to get the principial axis of inertia

      eigvals,eigvecs = LA.eig(all_mols[m].MIT)

# Transpose eigenvector matrix for vector transformation

      eigvecs=eigvecs.transpose()

# Now rotate all position vectors of the molecule along the prinicipial axis of inertia 
      all_mols[m].Rot_Molecule(eigvecs)
#
#   Now rotate the molecule to the desired orientation
#
      Emat=np.zeros((3,3))
      Emat=EulerMatrix(all_mols[m].prec,all_mols[m].nut,all_mols[m].rot)

      all_mols[m].Rot_Molecule(Emat)




#   Transform surface structure to xyz Angstrom coordinates
#   Transform them from direct coordinates by formula:
#   R = x1a1 + x2a2 + x3a3

   xyz_surf2 = np.zeros((natoms_surf,3))

   for i in range(natoms_surf): 
#   xyz_surf2[i][0]=(xyz_surf[i][0]*a_vec[0]+xyz_surf[i][1]*a_vec[1]+\
#                     xyz_surf[i][2]*a_vec[2])*surf_scale
#   xyz_surf2[i][1]=(xyz_surf[i][0]*b_vec[0]+xyz_surf[i][1]*b_vec[1]+\
#                     xyz_surf[i][2]*b_vec[2])*surf_scale
#   xyz_surf2[i][2]=(xyz_surf[i][0]*c_vec[0]+xyz_surf[i][1]*c_vec[1]+\
#                     xyz_surf[i][2]*c_vec[2])*surf_scale
      for j in range(3):
         xyz_surf2[i][j] = xyz_surf[i][j] * surf_scale 

# if surface rotation has been turned on, rotate it by 90 degrees
#   if rot_surf:
#      for i in range(natoms_surf):
#         x_act = xyz_surf2[i][0]
#         xyz_surf2[i][0] = xyz_surf[i][1]
#         xyz_surf2[i][1] = x_act

# Calculate COM of surface

   com_surf = COM(xyz_surf2,natoms_surf,names_surf,elements_dict,move=False)

   PrintMolecule(xyz_surf2,names_surf,natoms_surf,"surface","surface_test.xyz")
# Calculate highest z-value on surface
   surf_zmax = np.amax(xyz_surf2,axis=0)[2]
# Calculate lowest z-value on surface
   surf_zmin = np.amin(xyz_surf2,axis=0)[2]

   


# Allocate the lowest 40% of surface atoms (lowest 2 layers in the case of 5 layers)
# to fixed atoms 
   surf_fixed=[]
   for i in range(natoms_surf):
      if xyz_surf2[i][2] <= (surf_zmax-surf_zmin)*0.4:
         surf_fixed.append(True)
      else:
         surf_fixed.append(False)



# Translate x and y coordinates of molecule
#
#   Now combine xyz structures of surface and molecule
#   Depending on the number of molecules/species to be placed, predefine a number 
#   of rules for efficient placement..
#   Loop through the list of molecules to place them according to their 
#   preferred values 
#

#   Get total number of atoms
   natoms_tot = natoms_surf
#   Get total number of atoms belonging to adsorbated molecules 
   natoms_mols = 0
   for i in range(nspecs):
      natoms_tot = natoms_tot + all_mols[i].natoms 
      natoms_mols = natoms_mols + all_mols[i].natoms

   xyz_tot=np.zeros((natoms_tot,3))
   xyz_allmols=np.zeros((natoms_mols,3))
   names_tot=[]

   at_act = 0
   for j in range (nspecs):
      zmin_act = np.amin(all_mols[j].xyz,axis=0)[2]
      for i in range(all_mols[j].natoms):
         all_mols[j].xyz[i][0]=all_mols[j].xyz[i][0]+\
                 all_mols[j].xpos0
         all_mols[j].xyz[i][1]=all_mols[j].xyz[i][1]+\
                 all_mols[j].ypos0
         all_mols[j].xyz[i][2]=all_mols[j].xyz[i][2]+(surf_zmax-\
                 zmin_act+all_mols[j].surf_dist)
         xyz_allmols[at_act][0]=all_mols[j].xyz[i][0]
         xyz_allmols[at_act][1]=all_mols[j].xyz[i][1]
         xyz_allmols[at_act][2]=all_mols[j].xyz[i][2]
                
         at_act = at_act + 1 
#
#   If the cutout option is activated, remove all surface atoms that are 
#   too near to the adsorbates
#
   if cutout:
      dist_mat = distance_matrix(xyz_surf2,xyz_allmols,p=2)
      dist_mins = dist_mat.min(axis=1)

      min_mask=dist_mins < dist_remove
  
#
#    Filter the deleted atoms of the surface
#
      xyz_surf3=np.zeros((natoms_surf,3))
      elem_num_surf_new=[0]*nelem
      elem_num_del=[0]*nelem
      natoms_surf=0
      index=0
      index2=0

      names_surf=[]
      for i in range(nelem):
      # total number of atoms in the surface
         for j in range(elem_num_surf[i]):
             if not min_mask[index]:
                xyz_surf3[index2][0]=xyz_surf2[index][0]
                xyz_surf3[index2][1]=xyz_surf2[index][1]
                xyz_surf3[index2][2]=xyz_surf2[index][2]
                index2=index2+1
                elem_num_surf_new[i]=elem_num_surf_new[i]+1
                natoms_surf=natoms_surf+1
                names_surf.append(elements_surf[i])
             else: 
                elem_num_del[i]=elem_num_del[i] + 1 

             index=index+1
      
      print("Number of removed surface atoms due to adsorbate placement:")
      for i in range(nelem):
         if elem_num_del[i] > 0:
            print(" - " + elements_surf[i] + " : " + str(elem_num_del[i]))
#
#   If selective dynamics shall be used, choose the atoms to be moved 
#   based on their distance to the adsorbate(s)
#

   if select_dist:

      surf_fixed=[]
      count_true = 0

      dist_small=distance_matrix(xyz_allmols,xyz_surf2) < opt_dist/surf_scale
      for i in range(natoms_surf):
         counter = 0
         for j in range(natoms_mols):
            if dist_small[j][i]:
               counter = counter +1
         if counter > 0:
            surf_fixed.append(False)
            count_true = count_true +1
         else:
            surf_fixed.append(True)
      print (" The select_dist command was found.")
      print ("  " + str(count_true) + " atoms are activated for optimization.")


#
#    Indicate if a random placement shall be done
#
   if random:
      print (" The random command was found. Rotations will be chosen by chance,")
      print ("  the positions will be varied by +/- 0.5 Angstroms.")


 #  Combine structures of molecules and surface for global VASP input

   if cutout:
      for i in range(natoms_surf):
         names_tot.append(names_surf[i])
         for j in range(3):
            xyz_tot[i][j]=xyz_surf3[i][j]
      at_count = 0
   else:    
      for i in range(natoms_surf):
         names_tot.append(names_surf[i])
         for j in range(3):
            xyz_tot[i][j]=xyz_surf2[i][j]
      at_count = 0
   
   natoms_tot=natoms_surf

   for i in range(nspecs):
      natoms_tot=natoms_tot+all_mols[i].natoms 
      for j in range(all_mols[i].natoms): 
         names_tot.append(all_mols[i].names[j])
         surf_fixed.append(False)
         for k in range(3):
            xyz_tot[at_count+natoms_surf][k]=all_mols[i].xyz[j][k]
         at_count = at_count + 1
   PrintMolecule(xyz_tot,names_tot,natoms_tot,"surface+molecule","system_full.xyz")

#  Write POSCAR file for new combined system

   original_stdout=sys.stdout
   with open("POSCAR","w") as f:
      sys.stdout = f

   # rescale the z-component of the unit cell vector such that the desired 
   # total length in z-direction of the system is retained

      z_scale=z_length/(c_vec[2]*surf_scale)
     # c_vec[2]=c_vec[2]*z_scale

      print("New adsorbate system")
      print(1.0)
      print(str(a_vec[0]*surf_scale) + "  " + str(a_vec[1]*surf_scale) + "  " + str(a_vec[2]))
      print(str(b_vec[0]*surf_scale) + "  " + str(b_vec[1]*surf_scale) + "  " + str(b_vec[2]))
      print(str(c_vec[0]) + "  " + str(c_vec[1]) + "  " + str(c_vec[2]*surf_scale))

      new= False
      elements_all = []
      for i in range(natoms_tot):
         name_act = names_tot[i]
         new = True
         for j in range(i):
            if name_act == names_tot[j]:      
               new = False
         if new:
            elements_all.append(name_act) 
         new = False
     
      el_specs = len(elements_all)

      for i in range(el_specs):
         f.write(elements_all[i] + " ")
      f.write("\n")

      elem_num_all=[]
      for i in range(el_specs):
         elem_num_all.append(0)

      for i in range(natoms_tot):
         for j in range(el_specs):
            if names_tot[i] == elements_all[j]:
               elem_num_all[j] = elem_num_all[j]+1 

      
      for i in range(el_specs):
         f.write(str(elem_num_all[i]) + " ")
      f.write("\n")


      print("Selective Dynamics")
      print("Cartesian")


      for i in range(el_specs):
         for j in range(natoms_tot):
            if names_tot[j] == elements_all[i]:
               if surf_fixed[j]: 
                  print(str(xyz_tot[j][0]) + " " + str(xyz_tot[j][1]) + " " +
                       str(xyz_tot[j][2]) + select_false)
               else: 
                  print(str(xyz_tot[j][0]) + " " + str(xyz_tot[j][1]) + " " +
                       str(xyz_tot[j][2]) + select_true) 

sys.stdout=original_stdout

# Build the POTCAR file automatically by concatenating the single 
# element POTCAR files, remove old one if necessary
os.system("rm POTCAR")
os.system("touch POTCAR")
for elem in elements_all:
   os.system("cat " + potcar_path + elem + "/POTCAR >> POTCAR")
 
print('''
   build_adsorbate.py has succesfully finished!
   Positions written to POSCAR, PAW potentials written to POTCAR.
   Test xyz structure written to system_full.xyz.
   ''')
if zonly:
   print('''   The -zonly flag was present, only z-coordinates activated.
   ''')



