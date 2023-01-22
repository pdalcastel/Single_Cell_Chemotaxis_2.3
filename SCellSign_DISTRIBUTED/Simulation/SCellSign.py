from cc3d import CompuCellSetup
from math import pow
from math import sqrt

#global variables can be shared between Python modules
global random_seed
#random seed
random_seed = 1

global cellRad, cellVol
#cell radius
cellRad=10.0;
#cell volume as if it was a sphere
cellVol=4.19*cellRad*cellRad*cellRad

global J0, T, Lambda_F_actin
#J0 is the default energy value
J0=20.
#'Potts Temperature', also referred as 'average membrane fluctuations'
T=100. 
#Lambda_F_actin is the protrusion force coefficient. Negative sign makes the protrusion happen down gradient of FActin
Lambda_F_actin=-175.

global phiF, phiN, phiC 
#phiF is the Lamel volume fraction relative to the entire cell
phiF=0.05; 
#phiN for Nuc and phiC for Cyto
phiN=0.15; phiC=1.-phiN-phiF 

global Lx, Ly, Lz
#Lx Ly and Lz are the lattice sizes, they scale with the cell size and phiF to avoid percolation and compression by the lattice roof
Lx=int(4*(phiF*cellVol/3.14+((1-phiF)*cellVol*2/4.19)**(2/3))**(1/2))
Ly=Lx; Lz=int(2.1*cellRad)

global x1, y1
#initial cell position at the lattice center
x1=int(Lx/2.); y1=int(Ly/2.)

global VtF, VtN, VtC
#target volumes of each structure. The +0.5 creates aditional fluctuation around equilibrium value
VtF=phiF*cellVol+.5; VtN=phiN*cellVol+.5; VtC=phiC*cellVol+.5

global rho, delta, mu, chi
#rho normalizes the probability of Cyto conversion to Lamel
rho=0.5
#as delta increases, more conversions will happen, intensifying chemotaxis
delta=0.0
#as mu increases, more asymmetric is the Lamel creation, favoring chemotaxis. It saturates, so 10^6 and 10^7 will not behave different as 10^0 and 10^1 will.
mu=0.0
#chi is the tanh() multiplying factor
chi=1.0

global tSim, time_offset
#offset is the waiting time before data aquisition starts
time_offset = 0
#tSim is the simulation time
tSim = 100000+time_offset

#CC3D's "configureSimulation" function defines the model structure to be implemented in the simulation (Objects, Properties, Behaviors, Interactions, Dynamics, Initial Conditions and Boundary Conditions)
def configure_simulation():

    from cc3d.core.XMLUtils import ElementCC3D
    
    CompuCell3DElmnt=ElementCC3D("CompuCell3D",{"Revision":"20200821","Version":"4.2.3"})
    
    MetadataElmnt=CompuCell3DElmnt.ElementCC3D("Metadata")
    MetadataElmnt.ElementCC3D("NumberOfProcessors",{},"1")
    MetadataElmnt.ElementCC3D("DebugOutputFrequency",{},"1000")
    
    #specify basic Potts parameters: lattice size, boundary conditions, temperature, simulation time, copying neighbor order and random seed
    PottsElmnt=CompuCell3DElmnt.ElementCC3D("Potts")
    PottsElmnt.ElementCC3D("Dimensions",{"x":Lx,"y":Ly,"z":Lz})
    PottsElmnt.ElementCC3D("Steps",{},tSim)
    PottsElmnt.ElementCC3D("Temperature",{},T)
    PottsElmnt.ElementCC3D("RandomSeed",{},random_seed)
    PottsElmnt.ElementCC3D("NeighborOrder",{},"1")
    PottsElmnt.ElementCC3D("Boundary_x",{},"Periodic")
    PottsElmnt.ElementCC3D("Boundary_y",{},"Periodic")
    
    #initialize "cell types", NOT TO BE CONFUSED WITH THE CONCEPT OF A CELL PER SE. You rather refer to them as object types, as they will then compose a cell.
    PluginElmnt=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CellType"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"0","TypeName":"Medium"})
    PluginElmnt.ElementCC3D("CellType",{"Freeze":"","TypeId":"1","TypeName":"SUBS_A"})
    PluginElmnt.ElementCC3D("CellType",{"Freeze":"","TypeId":"2","TypeName":"SUBS_NA"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"3","TypeName":"CYTO"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"4","TypeName":"LAMEL"})
    PluginElmnt.ElementCC3D("CellType",{"TypeId":"5","TypeName":"NUC"})
    
    #Initialize target volumes and lambda volumes (inverse of compressibility) for each cell (object) type
    PluginElmnt_0=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Volume"})
    PluginElmnt_0.ElementCC3D("VolumeEnergyParameters",{"CellType":"CYTO","LambdaVolume":"10.0","TargetVolume":VtC})
    PluginElmnt_0.ElementCC3D("VolumeEnergyParameters",{"CellType":"LAMEL","LambdaVolume":"10.0","TargetVolume":VtF})
    PluginElmnt_0.ElementCC3D("VolumeEnergyParameters",{"CellType":"NUC","LambdaVolume":"10.0","TargetVolume":VtN})
    
    #Call important CC3D Plugins
    PluginElmnt_1=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
    PluginElmnt_2=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
    PluginElmnt_3=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"PixelTracker"})
    PluginElmnt_4=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"BoundaryPixelTracker"})
    PluginElmnt_4.ElementCC3D("NeighborOrder",{},"1")
    
    #Initialize contact energies between the compartments of DIFFERENT cells (per se)
    PluginElmnt_5=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Contact"})
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"Medium","Type2":"Medium"},1.)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"Medium","Type2":"SUBS_A"},J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"Medium","Type2":"SUBS_NA"},-J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"Medium","Type2":"CYTO"},J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"Medium","Type2":"LAMEL"},2.*J0/3)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"Medium","Type2":"NUC"},5.*J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"SUBS_A"},1.)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"SUBS_NA"},1.)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"CYTO"},J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"LAMEL"},J0/3)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"NUC"},5.*J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"SUBS_NA"},1.)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"CYTO"},J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"LAMEL"},2.*J0/3)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"NUC"},5.*J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"CYTO","Type2":"CYTO"},0.0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"CYTO","Type2":"LAMEL"},J0/2.)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"CYTO","Type2":"NUC"},J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"LAMEL","Type2":"LAMEL"},2.*J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"LAMEL","Type2":"NUC"},2.*J0)
    PluginElmnt_5.ElementCC3D("Energy",{"Type1":"NUC","Type2":"NUC"},0.0)
    PluginElmnt_5.ElementCC3D("NeighborOrder",{},"4")
    
    #Initialize contact energies between the compartments WITHIN a cell (per se)
    PluginElmnt_6=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"ContactInternal"})
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"SUBS_A"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"SUBS_NA"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"CYTO"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"LAMEL"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_A","Type2":"NUC"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"SUBS_NA"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"CYTO"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"LAMEL"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"SUBS_NA","Type2":"NUC"},"10.0")
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"CYTO","Type2":"CYTO"},0.)
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"CYTO","Type2":"LAMEL"},J0/2)
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"CYTO","Type2":"NUC"},J0)
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"LAMEL","Type2":"LAMEL"},0.0)
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"LAMEL","Type2":"NUC"},2.*J0)
    PluginElmnt_6.ElementCC3D("Energy",{"Type1":"NUC","Type2":"NUC"},0.)
    PluginElmnt_6.ElementCC3D("NeighborOrder",{},"4")
    
    #Initialize chemotaxis Plugin. THIS IS NOT THE CHEMOTAXIS PER SE, THIS IS FOR THE LAMEL PROTRUSION DYNAMICS
    PluginElmnt_7=CompuCell3DElmnt.ElementCC3D("Plugin",{"Name":"Chemotaxis"})
    ChemicalFieldElmnt=PluginElmnt_7.ElementCC3D("ChemicalField",{"Name":"FActin"})
    ChemicalFieldElmnt.ElementCC3D("ChemotaxisByType",{"ChemotactTowards":"Medium","Lambda":Lambda_F_actin,"Type":"LAMEL"})
    
    #Initialize diffusion solver for FActin field, but we do not use it, we update FActin values in the Steppables
    SteppableElmnt=CompuCell3DElmnt.ElementCC3D("Steppable",{"Type":"FlexibleDiffusionSolverFE"})
    DiffusionFieldElmnt=SteppableElmnt.ElementCC3D("DiffusionField",{"Name":"FActin"})
    DiffusionDataElmnt=DiffusionFieldElmnt.ElementCC3D("DiffusionData")
    DiffusionDataElmnt.ElementCC3D("FieldName",{},"FActin")
    DiffusionDataElmnt.ElementCC3D("DiffusionConstant",{},"0.0")
    DiffusionDataElmnt.ElementCC3D("DecayConstant",{},"0.0")
    DiffusionDataElmnt.ElementCC3D("ExtraTimesPerMCS",{},"0")
    
    '''FOR PIFF SAVING'''
    #SteppableElmnt_2=CompuCell3DElmnt.ElementCC3D("Steppable",{"Frequency":"100","Type":"PIFDumper"})
    # Periodically stores cell layout configuration in a piff format
    #SteppableElmnt_2.ElementCC3D("PIFName",{},"PiffLoaded")
    #SteppableElmnt_2.ElementCC3D("PIFFileExtension",{},"piff")

    CompuCellSetup.set_simulation_xml_description(CompuCell3DElmnt)

configure_simulation()

#next 3 lines execute the steppables' class that performs the model's dynamics
from SCellSignSteppables import Cell 
CellInstance=Cell(frequency=1,phiC=phiC,phiF=phiF,phiN=phiN,cellRad=cellRad,cellVol=cellVol,Lx=Lx,Ly=Ly,Lz=Lz,x1=x1,y1=y1,rho=rho,delta=delta,mu=mu,chi=chi,time_offset=time_offset)
CompuCellSetup.register_steppable(steppable=CellInstance)

'''uncomment next 3 lines to activate data output'''
#from SCellSignSteppables import Calc 
#CalcInstance=Calc(frequency=1,Lambda_F_actin=Lambda_F_actin,cellRad=cellRad,phiF=phiF,rho=rho,delta=delta,mu=mu,time_offset=time_offset,random_seed=random_seed)
#CompuCellSetup.register_steppable(steppable=CalcInstance)

CompuCellSetup.run()
