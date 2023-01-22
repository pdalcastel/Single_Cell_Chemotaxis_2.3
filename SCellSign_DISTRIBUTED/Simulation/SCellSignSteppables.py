#########################################################################################################################################
#    This code implements a simple model of a 3D crawling fibrobast-like cell with a nucleus (Nuc), cytoplasm (Cyto) 
# and lamellipodium (Lamel) that can respond to a chemical field biasing its migration direction (polarization).
#    The model for the cell uses a single chemical field (FActin) confined to the cell volume which represents the 
# region of actin polarization in the cell.
#    The chemical field is 1 in Lamel and 0 otherwise, and only in the Lamel pixels that make contact with Substrate.
#    The chemical field applies a repulsive chemotactic force of strength \lambda_F-actin (Eq. 5) to the boundary 
# between lamellipodium and medium only.
#    The cell starts migrating spontaneously after its shape symmetry is broken. 
#    The trajectory of the cell is calculated and recorded.
#    The code was originally written between 2014-2015 mainly by I. Fortuna and Gilberto L. Thomas (IFUFRGS, Brazil). 
#    It was then incrementally changed by Gabriel C. Perrone (IFUFRGS, Brazil) 2015-2016 aiming at reproducing chemotaxis.
#    Between 2018-2021 the code was uptaded and the chemotaxis model was finished by Pedro C. Dal Castel
#########################################################################################################################################

# Generic Python housekeeping steps -- loads utility functions
from math import pow
from math import sqrt
from math import tanh
from math import sin
from math import pi
from cc3d.core.PySteppables import *
import random
import os
import numpy

# Defines the Cell steppable class, which is called in SCellSign.py code
# In this class, we set the initial simulation state, solve chemotaxis dynmics, and update FActin field
class Cell(SteppableBasePy):

    def __init__(self,frequency,phiC,phiF,phiN,cellRad,cellVol,Lx,Ly,Lz,x1,y1,rho,delta,mu,chi,time_offset):

        SteppableBasePy.__init__(self,frequency)# Loads generic steppable class which executes at intervals specified by the CellMig3D.py file in mcs (here, once per mcs).
        # Imports parameters into "self" variables that can be used anywhere inside the class.
        self.phiF=phiF; self.phiC=phiC; self.phiN=phiN
        self.cellVol=cellVol; self.cellRad=cellRad
        self.Lx=Lx; self.Ly=Ly; self.Lz=Lz
        self.x1=x1; self.y1=y1
        self.rho=rho; self.delta=delta; self.mu=mu; self.chi=chi
        self.Saved_Lamel_Vol = []
        self.phiEST = self.phiF
        self.time_offset = time_offset
    
    def start(self):

        pt=CompuCell.Point3D()

        SUBS_Acell=self.potts.createCell() # Creates a cell (object)
        SUBS_Acell.type=self.SUBS_A # Gives a type to the created cell (object): Adherent Substrate
        SUBS_NAcell=self.potts.createCell() # Creates a cell (object)
        SUBS_NAcell.type=self.SUBS_NA # Gives a type to the created cell (object): Non Adherent Substrate
        CYTOcell=self.potts.createCell() # Creates a cell (object)
        CYTOcell.type=self.CYTO # Gives a type to the created cell (object): Cyto
        LAMELcell=self.potts.createCell() # Creates a cell (object)
        LAMELcell.type=self.LAMEL # Gives a type to the created cell (object): Lamel
        NUCcell=self.potts.createCell() # Creates a cell (object)
        NUCcell.type=self.NUC # Gives a type to the created cell (object): Nuc
        
        # The following routine draws the Cyto in a semisphere shape with the Nuc inside 
        # and Lamel surrounding the Cyto like an egg white
        Rad=int(self.cellRad)
        RadC=((1-self.phiF)*self.cellVol*2/4.19)**(1/3)
        RadF=(self.phiF*self.cellVol/3.14+((1-self.phiF)*self.cellVol*2/4.19)**(2/3))**(1/2)
        RadN=Rad*pow(self.phiN,1./3.)
        x1CM=self.x1
        y1CM=self.y1
        z1CM=Rad
        for x,y,z in self.every_pixel():
            pt.x=x; pt.y=y; pt.z=z
            if   ( pt.z==0 ):
                self.cellField.set(pt,SUBS_Acell)
            elif ( pt.z==self.dim.z-1 ):
                self.cellField.set(pt,SUBS_NAcell)
            elif ( pt.z==1 and sqrt((x-x1CM)*(x-x1CM)+(y-y1CM)*(y-y1CM)) <= RadF and sqrt((x-x1CM)*(x-x1CM)+(y-y1CM)*(y-y1CM)) > RadC ): 
                self.cellField.set(pt,LAMELcell)
            elif (sqrt((x-x1CM)*(x-x1CM)+(y-y1CM)*(y-y1CM)+(z-1)*(z-1)) < RadC): 
                if(sqrt((x-x1CM)*(x-x1CM)+(y-y1CM)*(y-y1CM)+(z-RadC/2)*(z-RadC/2)) <= RadN): 
                    self.cellField.set(pt,NUCcell) 
                else:
                    self.cellField.set(pt,CYTOcell) 

    def step(self,mcs):
        
        pt=CompuCell.Point3D() # Creates a vector variable to hold a position (x,y,z). 
        self.muef = self.mu # This does not do anything, see next line to understand how this would be useful
        # self.muef = self.mu*np.sin(mcs/10000) # This would model an temporally oscilating external field with period 10000 without the need to change the value of self.mu
        
        # The following routine updates the FActin every step: it sets FActin=1 for every Lamel pixels in z=1, and 0 elsewhere
        F = self.field.FActin
        F[:,:,1] = 0
        for cell in self.cell_list:
            if cell.type == 4:
                for pixel in self.get_cell_boundary_pixel_list(cell):
                    if pixel.pixel.z == 1:
                        F[pixel.pixel.x,pixel.pixel.y,pixel.pixel.z] = 1
        
        # Loop over Cyto type cells (objects)
        for cell in self.cellListByType(self.CYTO):  
            LAMELvol=0.0;   NUCvol=0.0  
            # Loop over Nuc type cells (objects) and get the current volume
            for NUCcell in self.cellListByType(self.NUC):
                NUCvol+=NUCcell.volume
            # Loop over Lamel type cells (objects) and get the current volume
            for LAMELcell in self.cellListByType(self.LAMEL):
                LAMELvol+=LAMELcell.volume 
            # Get total cell (per se) volume
            CELLvol = cell.targetVolume + LAMELvol + NUCvol 
            
            # Measure the average phiF over the previous 100 steps (Eq. 6)
            if mcs<100: 
                self.Saved_Lamel_Vol.append(LAMELvol/CELLvol)
            else:
                self.Saved_Lamel_Vol[mcs%100] = LAMELvol/CELLvol
                self.phiEST = sum(self.Saved_Lamel_Vol)/len(self.Saved_Lamel_Vol)

# CHEMOTAXIS      CHEMOTAXIS      CHEMOTAXIS      CHEMOTAXIS      CHEMOTAXIS      CHEMOTAXIS      CHEMOTAXIS 

            pList=[] 
            
            # Next 3 lines get all pixels that compose the cell(per se)'s surface and organize in the list pixelList
            cyto_pixelList = self.get_cell_boundary_pixel_list(cell)
            lamel_pixelList = self.get_cell_boundary_pixel_list(LAMELcell)
            pixelList = [[pixel.pixel.x,pixel.pixel.y,pixel.pixel.z] for pixel in cyto_pixelList]+[[pixel.pixel.x,pixel.pixel.y,pixel.pixel.z] for pixel in lamel_pixelList]
            
            chemosum=0.; chemosum2=0.; chemomed=0.; chemomed2=0.; normalization=0.; chemofactor=0.; chemomax=0.
            chemomin=0.; argchemomax=0.; argchemomin=0.; N=0.
            
            # Measure External Field concentration over cell(per se)'s base
            for pixel in pixelList: 
                if pixel[2]==1: # restrict to pixels at the base (z=1)
                    N+=1.0
                    pt.x=pixel[0]; pt.y=pixel[1]; pt.z=pixel[2]
                    which_cell = self.cellField.get(pt) # gets the cell to whom the pixel belong.
                    
                    # The field value is defined as the cell x position, but is corrected in case the cell is near a boundary
                    # This way, the cell always senses the same gradient throughout the simulation
                    if cell.xCOM > 2*self.dim.x/3 and pt.x < self.dim.x/3:
                        pt.x += self.dim.x
                    if cell.xCOM < self.dim.x/3 and pt.x > 2*self.dim.x/3:
                        pt.x -= self.dim.x
                    if (which_cell.type == 3): # gets all pixels in the base of Cyto
                        pList.append([pt.x,pt.y,pt.z]) # puts them in pList
                        
                    chemosum += pt.x
                    chemosum2 += pt.x**2
            chemomed = chemosum/N # average of External Field concentration
            chemomed2 = chemosum2/N 
            normalization = sqrt(chemomed2-chemomed*chemomed) # standard deviation of External Field concentration
            random.shuffle(pList)
            
            # Pixel conversions Cyto->Lamel
            for pixel in range(len(pList)): #Loop over all pixels from Cyto's base
                pt.x=pList[pixel][0]; pt.y=pList[pixel][1]; pt.z=pList[pixel][2]; 
                # Calculates the tanh() of the difference between the External Field concentration 
                # at the pixel and the average (Eq. 6)
                chemofactor = self.chi*tanh(self.muef*(pt.x-chemomed)/normalization)+1 
                # Weights the conversion probability by the tanh(), so that pixels sensing
                # higher External Field concentration will have higher chances of being converted to Lamel
                pLAMEL = self.rho*chemofactor 
                
                # Test if the current phiF is lower than the average over the last 100 mcs (Eq. 6)
                if (random.random()<pLAMEL and LAMELvol/CELLvol - self.delta <= self.phiEST):
                    
                    # Corrects x position in in case the cell is crossing a boundary
                    if pt.x >= self.dim.x:
                        pt.x -= self.dim.x
                    if pt.x < 0:
                        pt.x += self.dim.x
                    
                    # perform the conversion and updates the Lamel volume
                    self.cellField.set(pt,LAMELcell)
                    LAMELvol += 1

# This class is responsible for calculating and writing all cell (per se) compartments' centers of mass
class Calc(SteppableBasePy):
    
    def __init__(self,frequency,LambCHEM,cellRad,phiF,rho,delta,mu,time_offset,random_seed):
        SteppableBasePy.__init__(self,frequency) 
        self.LambCHEM = LambCHEM 
        self.cellRad = cellRad 
        self.phiF = phiF; self.rho = rho; self.delta = delta; self.mu = mu
        self.time_offset = time_offset
        self.random_seed = random_seed
  
    def start(self):
        
        out_dir_name = "SCellSign_output"
        print ("*************** AQUI ****************")
        print (out_dir_name)
        if not os.path.exists(out_dir_name): os.makedirs(out_dir_name)
        file_name = "_r"+str(self.cellRad)+"_f"+str(self.phiF)+"_lch"+str(self.LambCHEM)+"_mu"+str(self.mu)+"_d"+str(self.delta)+"_off"+str(self.time_offset)+"_rs"+str(self.random_seed)+"_Displacement.dat"
        self.output_path = str(Path(out_dir_name+"\\"+file_name))
        
        with open (self.output_path, 'w') as out_file:
            out_file.write("time   xposC        yposC         xposF          yposF          xposN          yposN          xposCN       yposCN         \n")
        out_file.close()
      
        self.posC=[.0,.0,.0]
        self.posF=[.0,.0,.0]
        self.posN=[.0,.0,.0]
        
        self.posrC=[.0,.0,.0]
        self.posrF=[.0,.0,.0]
        self.posrN=[.0,.0,.0]
        
        self.oposC=[.0,.0,.0]
        self.oposF=[.0,.0,.0]
        self.oposN=[.0,.0,.0]
        
        self.deslC=[.0,.0,.0]
        self.deslF=[.0,.0,.0]
        self.deslN=[.0,.0,.0]
        
        for cell in self.cellListByType(self.CYTO): 
            self.oposC=[cell.xCOM,cell.yCOM,cell.zCOM] 
            self.posrC=[cell.xCOM,cell.yCOM,cell.zCOM] 

        for cell in self.cellListByType(self.LAMEL):
            self.oposF=[cell.xCOM,cell.yCOM,cell.zCOM] 
            self.posrF=[cell.xCOM ,cell.yCOM,cell.zCOM] 

        for cell in self.cellListByType(self.NUC):
            self.oposN=[cell.xCOM,cell.yCOM,cell.zCOM] 
            self.posrN=[cell.xCOM,cell.yCOM,cell.zCOM] 

    def step(self,mcs):  
        
        if mcs >= self.time_offset: # Waiting before the calculation starts
            pt=CompuCell.Point3D() 

            # Calculate Cyto's center of mass and correct it in case it crossed a boundary
            for cell in self.cellListByType(self.CYTO): 

                self.posC[0]=cell.xCOM 
                self.posC[1]=cell.yCOM
                self.posC[2]=cell.zCOM
                self.deslC[0]=self.posC[0]-self.oposC[0]
                self.deslC[1]=self.posC[1]-self.oposC[1]
                
                if self.deslC[0] > 0.9*self.dim.x:
                    self.deslC[0] -= self.dim.x
                if self.deslC[0] < -0.9*self.dim.x:
                    self.deslC[0] += self.dim.x
                if self.deslC[1] > 0.9*self.dim.y:
                    self.deslC[1] -= self.dim.y
                if self.deslC[1] < -0.9*self.dim.y:
                    self.deslC[1] += self.dim.y

                self.posrC[0]=self.posrC[0]+self.deslC[0]
                self.posrC[1]=self.posrC[1]+self.deslC[1]
                self.oposC[0]=self.posC[0]
                self.oposC[1]=self.posC[1]
                self.oposC[2]=self.posC[2]
                
                self.volC = cell.volume 
            
            # Calculate Lamel's center of mass and correct it in case it crossed a boundary
            for cell in self.cellListByType(self.LAMEL): 

                self.posF[0]=cell.xCOM 
                self.posF[1]=cell.yCOM
                self.posF[2]=cell.zCOM
                self.deslF[0]=self.posF[0]-self.oposF[0] 
                self.deslF[1]=self.posF[1]-self.oposF[1] 

                if self.deslF[0] > 0.9*self.dim.x:
                    self.deslF[0] -= self.dim.x
                if self.deslF[0] < -0.9*self.dim.x:
                    self.deslF[0] += self.dim.x
                if self.deslF[1] > 0.9*self.dim.y:
                    self.deslF[1] -= self.dim.y
                if self.deslF[1] < -0.9*self.dim.y:
                    self.deslF[1] += self.dim.y

                self.posrF[0]=self.posrF[0]+self.deslF[0]
                self.posrF[1]=self.posrF[1]+self.deslF[1]
                self.oposF[0]=self.posF[0]
                self.oposF[1]=self.posF[1]
                self.oposF[2]=self.posF[2]
                
                self.volF = cell.volume 
            
            # Calculate Nuc's center of mass and correct it in case it crossed a boundary
            for cell in self.cellListByType(self.NUC): 

                self.posN[0]=cell.xCOM 
                self.posN[1]=cell.yCOM
                self.posN[2]=cell.zCOM
                self.deslN[0]=self.posN[0]-self.oposN[0] 
                self.deslN[1]=self.posN[1]-self.oposN[1] 

                if self.deslN[0] > 0.9*self.dim.x:
                    self.deslN[0] -= self.dim.x
                if self.deslN[0] < -0.9*self.dim.x:
                    self.deslN[0] += self.dim.x
                if self.deslN[1] > 0.9*self.dim.y:
                    self.deslN[1] -= self.dim.y
                if self.deslN[1] < -0.9*self.dim.y:
                    self.deslN[1] += self.dim.y

                self.posrN[0]=self.posrN[0]+self.deslN[0]
                self.posrN[1]=self.posrN[1]+self.deslN[1]
                self.oposN[0]=self.posN[0]
                self.oposN[1]=self.posN[1]
                self.oposN[2]=self.posN[2]
                
                self.volN = cell.volume 
            
            # Normalize positions by cell radius
            xC=self.posrC[0]/self.cellRad
            yC=self.posrC[1]/self.cellRad
            zC=self.posC[2]/self.cellRad
            xF=self.posrF[0]/self.cellRad
            yF=self.posrF[1]/self.cellRad
            zF=self.posF[2]/self.cellRad
            xN=self.posrN[0]/self.cellRad
            yN=self.posrN[1]/self.cellRad
            zN=self.posN[2]/self.cellRad
            
            # Calculate position of Cyto and Nuc together
            xCN=(self.volN*xN+self.volC*xC)/(self.volN+self.volC)
            yCN=(self.volN*yN+self.volC*yC)/(self.volN+self.volC)
            
            FCV = 1.*self.volF/(self.volC+self.volN+self.volF) 
            
            # Write output into a file
            with open (self.output_path, 'a') as out_file:
                out_file.write("%d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n" % (mcs-self.time_offset ,xC,yC,xF,yF,xN,yN,xCN,yCN))
            out_file.close()
