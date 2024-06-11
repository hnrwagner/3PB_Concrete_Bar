# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 12:41:49 2022

@author: Ronald Wagner
"""

from abaqus import *
from abaqusConstants import *
import regionToolset
import __main__
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import odbAccess
from operator import add


import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# functions

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def Create_3D_Beam(model,part,length,height,thickness):
    s = mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(-length/2.0, height/2.0), point2=(length/2.0, -height/2.0))
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models[model].parts[part]
    p.BaseSolidExtrude(sketch=s, depth=thickness)
    s.unsetPrimaryObject()
    del mdb.models[model].sketches['__profile__']

#------------------------------------------------------------------------------    

def Create_Datum_Plane_by_Principal(type_plane,part,model,offset_plane):  
    p = mdb.models[model].parts[part]
    myPlane = p.DatumPlaneByPrincipalPlane(principalPlane=type_plane, offset=offset_plane)
    myID = myPlane.id
    return myID    

#------------------------------------------------------------------------------

def Create_Partion_by_Plane(model,part,id_plane):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    d = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d[id_plane], cells=c)
    
#------------------------------------------------------------------------------

def Create_Set_All_Cells(model,part,set_name):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    p.Set(cells=c, name=set_name) 

#------------------------------------------------------------------------------

def Create_Set_Face(x,y,z,model,part,set_name):
    face = ()
    p = mdb.models[model].parts[part]
    f = p.faces
    myFace = f.findAt((x,y,z),)
    face = face + (f[myFace.index:myFace.index+1],)
    p.Set(faces=face, name=set_name)
    return myFace

#------------------------------------------------------------------------------

def Create_Assembly(model,part,instance,x,y,z):
    a = mdb.models[model].rootAssembly
    p = mdb.models[model].parts[part]
    a.Instance(name=instance, part=p, dependent=ON)
    p = a.instances[instance]
    p.translate(vector=(x,y,z))
    
#------------------------------------------------------------------------------

def Create_Reference_Point(x,y,z,model,setname):
    a = mdb.models[model].rootAssembly
    myRP = a.ReferencePoint(point=(x, y, z))
    r = a.referencePoints
    myRP_Position = r.findAt((x, y, z),)    
    refPoints1=(myRP_Position, )
    a.Set(referencePoints=refPoints1, name=setname)
    return myRP,myRP_Position

#------------------------------------------------------------------------------

def Create_Interaction_Coupling(model,instance,rp_name,rp_face_name,constraint_name):
    a = mdb.models[model].rootAssembly
    region1=a.sets[rp_name]
    region2=a.instances[instance].sets[rp_face_name]
    mdb.models[model].Coupling(name=constraint_name, controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=STRUCTURAL, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

#------------------------------------------------------------------------------

def Create_Analysis_Step(model,step_name,pre_step_name,Initial_inc,Max_inc,Min_inc,Inc_Number,NL_ON_OFF):
    a = mdb.models[model].StaticStep(name=step_name, previous=pre_step_name, initialInc=Initial_inc, maxInc=Max_inc, minInc=Min_inc, nlgeom=ON)
    a = mdb.models[model].steps[step_name].setValues(maxNumInc=Inc_Number)
    a = mdb.models[model].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CSTRESS', 'CDISP', 'DAMAGEC', 'DAMAGET','STATUS'))
    a = mdb.models[model].steps[step_name].setValues(stabilizationMagnitude=1E-009, stabilizationMethod=DAMPING_FACTOR, continueDampingFactors=False, adaptiveDampingRatio=None)
 
#------------------------------------------------------------------------------

def Create_Gravity_Load(model,load_name,step_name,load):
    mdb.models[model].Gravity(name=load_name, createStepName=step_name, comp2=-load, distributionType=UNIFORM, field='')

#------------------------------------------------------------------------------

def Create_BC(model,rp_name,BC_name,step_name,u,v,w,ur,vr,wr):
    a = mdb.models[model].rootAssembly
    region = a.sets[rp_name]
    mdb.models[model].DisplacementBC(name=BC_name, createStepName=step_name, region=region, u1=u, u2=v, u3=w, ur1=ur, ur2=vr, ur3=wr, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

#------------------------------------------------------------------------------

def Create_BC_2(model,instance_name,set_name,BC_name,step_name,u,v,w,ur,vr,wr):
    a = mdb.models[model].rootAssembly
    region = a.instances[instance_name].sets[set_name]
    mdb.models[model].DisplacementBC(name=BC_name, createStepName=step_name, region=region, u1=u, u2=v, u3=w, ur1=ur, ur2=vr, ur3=wr, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

#------------------------------------------------------------------------------

def Create_Mesh(model,part,mesh_size):
    p = mdb.models[model].parts[part]
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

#------------------------------------------------------------------------------

def Create_Material_and_Assign(model,part,material_name,E,Nu,Rho,section_name,set_name):
    p = mdb.models[model].parts[part]
    mdb.models[model].Material(name=material_name)
    mdb.models[model].materials[material_name].Elastic(table=((E, Nu), ))
    mdb.models[model].materials[material_name].Density(table=((Rho, ), ))
    mdb.models[model].HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)
    p = mdb.models[model].parts[part]
    region = p.sets[set_name]
    p = mdb.models[model].parts[part]
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    
#------------------------------------------------------------------------------

def Create_Job(model,job_name, cpu):
    a = mdb.models[model].rootAssembly
    mdb.Job(name=job_name, model=model, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=cpu, numDomains=cpu, numGPUs=0)

#------------------------------------------------------------------------------

def SubmitJob(job_name):
    mdb.jobs[job_name].submit(consistencyChecking=OFF)
    mdb.jobs[job_name].waitForCompletion()
#------------------------------------------------------------------------------    

def Open_ODB_and_Write_NodeSet_data_to_text(model,step_name,variable_name,set_name,Variable_component):
    # open ODB file - ABAQUS Result file
    odb = session.openOdb(str(model)+'.odb')
    
    # list for the VARIABLE you want to evaluate
    Variable_v = []
    
    # analysis step for your VARIABLE
    lastStep=odb.steps[step_name]
    
    #loop over all increments of the analysis step and save VARIABLE information from each increment
    for x in range(len(lastStep.frames)):
        lastFrame = lastStep.frames[x]
        Variable = lastFrame.fieldOutputs[variable_name]
        center = odb.rootAssembly.nodeSets[set_name]
        centerRForce = Variable.getSubset(region=center)
       
        # loop over the VARIABLE and save component (x,y,z - 0,1,2) to list
        for i in centerRForce.values:
            Variable_vr = [i.data[Variable_component]]
            Variable_v = Variable_v + Variable_vr  
            
    # write VARIABLE - component to text file
    
    np.savetxt(str(set_name)+'_'+str(variable_name)+'_'+str(myString)+'.txt',Variable_v)
    return Variable_v

#------------------------------------------------------------------------------  

def Create_Sum_ABS_List(data_1,data_2,factor):
    data = [(x + y)/(2.0*factor) for (x, y) in zip(data_1, data_2)] 
    data =  [abs(ele) for ele in data]
    max_data = max(data)
    return max_data,data

#------------------------------------------------------------------------------  

def Create_Plot(data_1,data_2,max_data_1,max_data_2):
    fig, ax = plt.subplots()
    ax.plot(data_1, data_2, 'g--', label='Concrete')
    legend = ax.legend(loc='lower right', shadow=True)
    plt.ylabel('Force [N]')
    plt.xlabel('Displacement [mm]')
    plt.grid(True)
    plt.axis([0.0, max_data_1, 0.0, max_data_2])
    plt.savefig('Load_Displacement_Curve_Concrete_3PB.png')

#------------------------------------------------------------------------------

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

def data_to_tuples_of_tuples(data1 , data2):
    list_of_list = [[i, j] for i,j in zip(data1 , data2)]
    tuples_of_tuples = tuple(map(tuple, list_of_list))
    return tuples_of_tuples

#---------------------------------------------------------------------------------------------------------------------------------

def Create_CDP(model,part,material_name,dil_angle,eccentricity,fb0_fc0,K,visc_para,CDP_Comp_Yield_Inelastic,CDP_Tension_Yield_Crack,CDP_Comp_Damage_Inelastic,CDP_Tension_Damage_Crack):
    p = mdb.models[model].parts[part]
    myCDPMaterial = mdb.models[model].materials[material_name]
    myCDPMaterial.ConcreteDamagedPlasticity(table=((dil_angle,eccentricity,fb0_fc0,K,visc_para), ))
    myCDPMaterial.concreteDamagedPlasticity.ConcreteCompressionHardening(table=CDP_Comp_Yield_Inelastic)
    myCDPMaterial.concreteDamagedPlasticity.ConcreteTensionStiffening(table=CDP_Tension_Yield_Crack)
    myCDPMaterial.concreteDamagedPlasticity.ConcreteCompressionDamage(table=CDP_Comp_Damage_Inelastic)
    myCDPMaterial.concreteDamagedPlasticity.ConcreteTensionDamage(table=CDP_Tension_Damage_Crack)

#---------------------------------------------------------------------------------------------------------------------------------

def Create_CS(model,part,radius,length):
    s = mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(radius, 0.0))
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    p = mdb.models[model].parts[part]
    p.BaseSolidExtrude(sketch=s, depth=length)
    s.unsetPrimaryObject()
    p = mdb.models[model].parts[part]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models[model].sketches['__profile__']
    
    # Convert Cell to Shell
    
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    p.RemoveCells(cellList = c)
    
    # Create Reference Point in middle of Support
    
    v, e, d, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=p.InterestingPoint(edge=e[0], rule=CENTER))
    
    # Create Set for Boundary Conditions
    
    c = p.cells[:]
    region = p.Set(cells=c, name='Set-CS')

#---------------------------------------------------------------------------------------------------------------------------------

def Material_CS(model,part):
    mdb.models[model].Material(name='Steel')
    mdb.models[model].materials['Steel'].Elastic(table=((210000.0, 0.3), ))
    mdb.models[model].materials['Steel'].Density(table=((7.85E-09, ), ))
    mdb.models[model].HomogeneousSolidSection(name='Section-CS', material='Steel', thickness=None)
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    region = p.Set(cells=c, name='Set-CS')
    #p = mdb.models[model].parts[part]
    #p.SectionAssignment(region=region, sectionName='Section-CS', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

#---------------------------------------------------------------------------------------------------------------------------------

def Assembly_CS(model,part,x,y,z,nr):
    a = mdb.models[model].rootAssembly
    p = mdb.models[model].parts[part]
    a.Instance(name='Cylinder_Support-'+str(nr), part=p, dependent=ON)
    p1 = a.instances['Cylinder_Support-'+str(nr)]
    p1.translate(vector=(x, y, z))
    
#---------------------------------------------------------------------------------------------------------------------------------

def Create_Contact(model):
    a = mdb.models[model].rootAssembly
    mdb.models[model].ContactProperty('IntProp-1')
    mdb.models[model].interactionProperties['IntProp-1'].TangentialBehavior(formulation=FRICTIONLESS)
    mdb.models[model].interactionProperties['IntProp-1'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    mdb.models[model].ContactStd(name='Int-2', createStepName='Initial')
    mdb.models[model].interactions['Int-2'].includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
    mdb.models[model].interactions['Int-2'].contactPropertyAssignments.appendInStep(stepName='Initial', assignments=((GLOBAL, SELF, 'IntProp-1'), ))

#---------------------------------------------------------------------------------------------------------------------------------

def Create_Rigid_Body(model,csnr,rpnr,cnr):
    a = mdb.models[model].rootAssembly
    region4=a.instances['Cylinder_Support-'+str(csnr)].sets['Set-CS']
    a = mdb.models[model].rootAssembly
    region1=a.sets['RP-'+str(rpnr)]
    mdb.models[model].RigidBody(name='Constraint-'+str(cnr), refPointRegion=region1, tieRegion=region4)
    

#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

# variables

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------    

myPart = "Frame"
myString = "Concrete_Frame_3PB_with_Damage-w5"
mdb.Model(name=myString)
myMaterialName = "Concrete"
myJobName= "3PB_Concrete_Beam_with_Damage-w5"


myLength = 750.0
myHeight = 150.0
myThickness = 150.0

myLoad = 5.0


myBCLength = 50.0


myE = 35500
myNu = 0.2
myRho = 2.5E-09

myE_Compression = 34007.32018
myY_Compression = 42
myE_Tension = 37000
myY_Tension = 2.0

myDila_Angle = 35
myEccentricity = 0.1
myFb0_Fc0 = 1.16
myK = 0.67
myVisco_Para = 1E-05

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# CDP - Start

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# Open Txt file with Stress-Strain data - Compression and Tension - Concrete

myCDP_Compression_Strain_list, myCDP_Compression_Stress_list = Open_txt_to_List("Stress_Strain_Compression.txt")
myCDP_Tension_Strain_list, myCDP_Tension_Stress_list = Open_txt_to_List("Stress_Strain_Tension.txt")

# Determine Length of Stress-Strain list

myIndex = len(myCDP_Compression_Stress_list)
myIndex2 = len(myCDP_Tension_Stress_list)

# Determine Damage Parameter, Inelastic Strain, Cracking Strain, Yield Stress

myDamageParameter_Compression, myInelastic_Compression_Strain, myCDP_Compression_Stress_list = Stress_Strain_to_ABAQUS_CDP(myIndex,myCDP_Compression_Strain_list,myCDP_Compression_Stress_list,myE_Compression,myY_Compression)
myDamageParameter_Tension, myCracking_Tension_Strain, myCDP_Tension_Stress_list = Stress_Strain_to_ABAQUS_CDP(myIndex2,myCDP_Tension_Strain_list,myCDP_Tension_Stress_list,myE_Tension,myY_Tension)

# Create data for Strain (Inelastic & Cracking) and Stress of Compression and Tension

myCDP_Comp_Yield_Inelastic = data_to_tuples_of_tuples(myCDP_Compression_Stress_list, myInelastic_Compression_Strain)
myCDP_Tension_Yield_Crack = data_to_tuples_of_tuples(myCDP_Tension_Stress_list, myCracking_Tension_Strain)

# Create data for Damage Parameter of Compression and Tension

myCDP_Comp_Damage_Inelastic = data_to_tuples_of_tuples(myDamageParameter_Compression, myInelastic_Compression_Strain)
myCDP_Tension_Damage_Crack = data_to_tuples_of_tuples(myDamageParameter_Tension, myCracking_Tension_Strain)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# CDP - End

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


# Create 3D Beam

Create_3D_Beam(myString, myPart, myLength, myHeight, myThickness)

# Create Planes for Partitions

myID_0 = Create_Datum_Plane_by_Principal(YZPLANE,myPart,myString,0.0+myBCLength/2.0)
myID_1 = Create_Datum_Plane_by_Principal(YZPLANE,myPart,myString,0.0-myBCLength/2.0)

myID_2 = Create_Datum_Plane_by_Principal(YZPLANE,myPart,myString,-225.0+myBCLength/2.0)
myID_3 = Create_Datum_Plane_by_Principal(YZPLANE,myPart,myString,-225.0-myBCLength/2.0)

myID_4 = Create_Datum_Plane_by_Principal(YZPLANE,myPart,myString,225.0+myBCLength/2.0)
myID_5 = Create_Datum_Plane_by_Principal(YZPLANE,myPart,myString,225.0-myBCLength/2.0)

# Create Partitions

Create_Partion_by_Plane(myString,myPart,myID_0)
Create_Partion_by_Plane(myString,myPart,myID_1)
Create_Partion_by_Plane(myString,myPart,myID_2)
Create_Partion_by_Plane(myString,myPart,myID_3)
Create_Partion_by_Plane(myString,myPart,myID_4)
Create_Partion_by_Plane(myString,myPart,myID_5)

# Create Set - Cell for whole model

Create_Set_All_Cells(myString,myPart,"Beam_3D")

# Create Set - Face RP-1 & RP-2

Create_Set_Face(0.0,myHeight/2.0,myThickness/2.0,myString,myPart,"Face_RP_1")

# Creaete Set - Face RP-3 & RP-4

Create_Set_Face(-225.0,-myHeight/2.0,myThickness/2.0,myString,myPart,'Face_RP_3')
Create_Set_Face(225.0,-myHeight/2.0,myThickness/2.0,myString,myPart,'Face_RP_4')

# Create Assembly

Create_Assembly(myString,myPart,"Beam-1",0,0,0)

# Create Analysis Step

Create_Analysis_Step(myString,"Loading","Initial",0.001,0.01,1E-15,300,ON)    


# Create Material Definition

Create_Material_and_Assign(myString,myPart,myMaterialName,myE,myNu,myRho,"Concrete_Beam_3D_Section","Beam_3D")

# Create CDP Materail Definition

Create_CDP(myString,myPart,"Concrete",myDila_Angle,myEccentricity,myFb0_Fc0,myK,myVisco_Para,myCDP_Comp_Yield_Inelastic,myCDP_Tension_Yield_Crack,myCDP_Comp_Damage_Inelastic,myCDP_Tension_Damage_Crack)

# Creaete Set - RP RP-1

myRP1,myRP_Position1 = Create_Reference_Point(0.0,myHeight/2.0,myThickness/2.0,myString,'RP-1')

# Create Coupling Interaction between RP-1 & RP-2 and Corresponding Faces

for ic in range(1,2,1):
    Create_Interaction_Coupling(myString,"Beam-1","RP-"+str(ic),"Face_RP_"+str(ic),"RP-"+str(ic)+"_to_Face")

# Create Boundary Conditions

# Create Loading
Create_BC(myString,"RP-1","Loading","Loading",SET,-myLoad,SET,SET,SET,SET)

Create_BC_2(myString,"Beam-1",'Face_RP_3',"Support_1","Initial",SET,SET,SET,SET,SET,SET)
Create_BC_2(myString,"Beam-1",'Face_RP_4',"Support_2","Initial",SET,SET,SET,SET,SET,SET)

# Create Mesh

Create_Mesh(myString,myPart,10.0)

# Create Job

Create_Job(myString, myJobName, 4)

# Submit Job

#SubmitJob(myJobName)

# Create Txt File with Load-Displacement Data

RF_data_1 = Open_ODB_and_Write_NodeSet_data_to_text(myJobName,'Loading','RF','RP-1',1)
U_data_1 = Open_ODB_and_Write_NodeSet_data_to_text(myJobName,'Loading','U','RP-1',1)

# # Create Figure of Load-Displacement Data

Len_RF, RF_data = Create_Sum_ABS_List(RF_data_1,RF_data_1,1)
Len_U, U_data = Create_Sum_ABS_List(U_data_1,U_data_1,1)

Create_Plot(U_data,RF_data,Len_U,Len_RF)


