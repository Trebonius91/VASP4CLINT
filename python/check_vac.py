#!/usr/bin/env python3.6


in_file = open('check_vac.in','r')
in_list = in_file.readlines()
in_file.close()

elements = in_list[0].split()
vac = in_list[1].split()
thresholds = in_list[2].split()

#print(elements)
#print(vac)
#print(thresholds)


contcar = open('CONTCAR', 'r')
contcar_list = contcar.readlines()
contcar.close()


element_list = contcar_list[5].split()
number_list = contcar_list[6].split()

#print(number_list)
#print(element_list)

ele_num = []
for x in range(len(element_list)):
    dummy = [element_list[x],int(number_list[x])]
    ele_num.append(dummy)
    #print(ele_num)
    

coord  = contcar_list[7].strip()
#print(coord)


lower_vac = float(vac[0]) + float(thresholds[0])
upper_vac = float(vac[1]) - float(thresholds[1])

#print(lower_vac)
#print(upper_vac)

Ntot = 0
for x in range(len(ele_num)):
    Ntot = Ntot + ele_num[x][1]
    #print(Ntot)

gas = False
for e in range(len(element_list)):
    if element_list[e] in elements:

        lower = 8
        for i in range(e):
            lower = lower + ele_num[i][1]

        #print(lower)
        #print(contcar_list[lower])
        #print(contcar_list[lower+ele_num[e][1]-1])

        for r in range(lower,lower+ele_num[e][1],1):

            if coord == "Direct":
                rz = float(contcar_list[r].split()[2]) * float(vac[1])       
            elif coord == "Cartesian":
                rz = float(contcar_list[r].split()[2])
            
            if rz < 0.0:
                rz = float(vac[1]) + rz


            #print(r , rz)

            if rz > lower_vac and rz < upper_vac:
                gas = True 
                #print(r-7, "danger")

if not gas:
    cont_scratch = open("CONTCAR_scratch",'w')
    for x in range(Ntot+8):
        cont_scratch.write(contcar_list[x])
    cont_scratch.close()
elif gas:
    kill = open("gas_kill",'w')
    kill.write('Gasphase-Atom detected')
    kill.close()
#print(gas)

