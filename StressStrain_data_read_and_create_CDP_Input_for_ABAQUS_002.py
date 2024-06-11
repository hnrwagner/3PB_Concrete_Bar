# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:17:11 2022

@author: Besitzer
"""


import numpy as np
import matplotlib.pyplot as plt

def index_of_max_in_list(list_data, max_value):
    for i in list_data:
        if i >max_value: break
    return list_data.index(i)

#---------------------------------------------------------------------------------------------------------------------------------

def Stress_Strain_to_ABAQUS_CDP(index_list,CDP_Load_Strain_list,CDP_Load_Stress_list,E_Load,Y_Load):
    Inelastic_Compression_Strain = []
    
    for ic in range(0,index_list,1):
        Inelastic_Compression_Strain = CDP_Load_Strain_list-CDP_Load_Stress_list/E_Load
            
        
    Inelastic_Compression_Strain = Inelastic_Compression_Strain.tolist()
    CDP_Load_Stress_list = CDP_Load_Stress_list.tolist()
    
    DamageParameter_Compression = []
    for Stress in CDP_Load_Stress_list:
        DamageParameter_Compression.append(1.0-(Stress/Y_Load))
       

    MinIndex_2 = index_of_max_in_list(CDP_Load_Stress_list,Y_Load)

    
    del DamageParameter_Compression[0:MinIndex_2]
    del Inelastic_Compression_Strain[0:MinIndex_2]
    del CDP_Load_Stress_list[0:MinIndex_2]
    
    MinIndex= DamageParameter_Compression.index(min([i for i in DamageParameter_Compression if i > 0]))
    
    DamageParameter_Compression = np.array(DamageParameter_Compression)
    DamageParameter_Compression[0:MinIndex] = 0
    DamageParameter_Compression = DamageParameter_Compression.tolist()
    Inelastic_Compression_Strain[0] = 0
    return DamageParameter_Compression, Inelastic_Compression_Strain, CDP_Load_Stress_list
    
#---------------------------------------------------------------------------------------------------------------------------------

def Open_txt_to_List(txt_name):
    CDP_Load_Stress_Strain_txt = np.loadtxt(txt_name)[:]
    CDP_Load_Strain_list = []
    CDP_Load_Strain_list = CDP_Load_Stress_Strain_txt[:,0]

    CDP_Load_Stress_list = []
    CDP_Load_Stress_list = CDP_Load_Stress_Strain_txt[:,1]
    return  CDP_Load_Strain_list,CDP_Load_Stress_list

#---------------------------------------------------------------------------------------------------------------------------------

def Create_XY_Plot(X,Y,name,X_name,Y_name,X_max,Y_max):    
    fig, ax = plt.subplots()
    ax.plot(X, Y, 'g--', label=name)
    legend = ax.legend(loc='lower right', shadow=True)
    plt.ylabel(Y_name)
    plt.xlabel(X_name)
    plt.grid(True)
    plt.axis([0.0, X_max, 0.0, Y_max])

#---------------------------------------------------------------------------------------------------------------------------------


myE_Compression = 34007.32018
myY_Compression = 42
myE_Tension = 37000
myY_Tension = 2.0


myCDP_Compression_Strain_list, myCDP_Compression_Stress_list = Open_txt_to_List("Stress_Strain_Compression.txt")
myCDP_Tension_Strain_list, myCDP_Tension_Stress_list = Open_txt_to_List("Stress_Strain_Tension.txt")

myIndex = len(myCDP_Compression_Stress_list)
myIndex2 = len(myCDP_Tension_Stress_list)

myDamageParameter_Compression, myInelastic_Compression_Strain, myCDP_Compression_Stress_list = Stress_Strain_to_ABAQUS_CDP(myIndex,myCDP_Compression_Strain_list,myCDP_Compression_Stress_list,myE_Compression,myY_Compression)
myDamageParameter_Tension, myCracking_Tension_Strain, myCDP_Tension_Stress_list = Stress_Strain_to_ABAQUS_CDP(myIndex2,myCDP_Tension_Strain_list,myCDP_Tension_Stress_list,myE_Tension,myY_Tension)


Create_XY_Plot(myInelastic_Compression_Strain,myCDP_Compression_Stress_list,"concrete","Inelastic Strain","Stress",max(myInelastic_Compression_Strain),max(myCDP_Compression_Stress_list))
Create_XY_Plot(myInelastic_Compression_Strain,myDamageParameter_Compression,"concrete","Inelastic Strain","Damage Parameter",max(myInelastic_Compression_Strain),1.0)

Create_XY_Plot(myCracking_Tension_Strain,myCDP_Tension_Stress_list,"concrete","Cracking Strain","Stress",max(myCracking_Tension_Strain),max(myCDP_Tension_Stress_list))
Create_XY_Plot(myCracking_Tension_Strain,myDamageParameter_Tension,"concrete","Cracking Strain","Damage Parameter",max(myCracking_Tension_Strain),1.0)

np.savetxt("Damage.txt",myDamageParameter_Compression,fmt='%1.4e')
np.savetxt("Strain.txt",myInelastic_Compression_Strain,fmt='%1.4e')



