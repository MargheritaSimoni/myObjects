
import NCrystal as NC
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

"""ELEMENT CLASS"""
class Element:
  def __init__(self, name, mass, sigmaBound):
        self.name = name
        self.mass = mass
        self.sigmaBound = sigmaBound
  def sigmaFree(self):
    neutronMass=1.008664
    sigmaFree=self.sigmaBound/((1+(neutronMass/self.mass))**2)
    return sigmaFree


"""MATERIAL CLASS"""
class Material:
    def __init__(self, name, formula, density, ncmat="AFGA", functionalGroups=None, temperature=293.6, elements_dictionary=False):
        self.name = name
        self.formula = formula
        self.density = density
        self.temperature = temperature
        self.ncmat = ncmat
        self.functionalGroups = functionalGroups
        self.elements_dictionary=elements_dictionary
        if elements_dictionary==False:
          raise ValueError(f"Elements_dictionary needs to be given as an argument")
        """build ncmat file"""
        if ncmat=="AFGA":
          if functionalGroups==None:
            print("ERROR: When using AFGA model, functional groups must be declared!")
          import subprocess
          command = ['ncrystal_hfg2ncmat','--formula', formula,'--spec', functionalGroups,'--density', str(density),'--title', name,'-o', f'{name}.ncmat',"--force"]
          creating_ncmat = subprocess.run(command, capture_output=True, text=True)
          self.creating_ncmat= creating_ncmat
          T=str(temperature)+"K"
          fn = f'{name}.ncmat;temp='+T
        else:
          #print("WARNING: Density and temperature may be defined differently inside NCrystal .ncmat file")
          fn=ncmat
        self.fn=fn
        """Define material stoichiometry from formula"""
         #______________________________________________#
        def define_Stoichiometry(formula):
          import re
          #Regular expression pattern to match elements and their quantities
          pattern = r'([A-Z][a-z]?)(\d*)'
          # Find all matches in the formula string
          matches = re.findall(pattern, formula)
          elements = []
          stoichiometry = []
          for (element, pedice) in matches:
              elements.append(element)
              # If the quantity is an empty string, it means there is only one atom of that element
              if pedice == '':
                  stoichiometry.append(1)
              else:
                  stoichiometry.append(int(pedice))
          return elements, stoichiometry
        elements, stoichiometry = define_Stoichiometry(formula)
        self.elements= elements
        self.stoichiometry =stoichiometry
        for el in self.elements:
          if el not in self.elements_dictionary: raise ValueError(f"Element {el} not found in elements_dictionary")

        """SIGMAFREE OF MATERIAL"""
        SigmaFree=0.
        for el, stoy in zip(self.elements, self.stoichiometry):
          SigmaFree+=elements_dictionary[el].sigmaFree()*stoy
        self.sigmaFree=SigmaFree
        """MASS OF MATERIAL"""
        mass=0.
        for el, stoy in zip(self.elements, self.stoichiometry):
            mass+=elements_dictionary[el].mass*stoy
        self.mass=mass


    """END OF CONSTRUCTOR"""
    #______________________________#
    def __str__(self):
        if self.ncmat=="AFGA":
          return f"\n_____________________\nMATERIAL PROPERTIES:\n_____________________\n name={self.name} \n formula={self.formula} \n density={self.density} \n functionalGroups={self.functionalGroups} \n temperature={self.temperature}\n\n This material was created using the Average functional group approximation (AFGA):\n {self.creating_ncmat.stdout}\n {self.creating_ncmat.stderr}\n"
        else:
          return f"\n_____________________\nMATERIAL PROPERTIES:\n_____________________\n name={self.name}\n formula={self.formula}\n density={self.density}\n ncmat={self.ncmat}\n temperature={self.temperature}\n\n This material was created using NCrystal libraries\n WARNING: Density and temperature may be defined differently inside NCrystal .ncmat file\n"


    """MICROSCOPIC CROSS SECTION OF MATERIAL"""
     #________________________________________#
    def crossSection(self, E, absorption=True, formulaunit=True, plot=False, help=False):
      pc = NC.createScatter(self.fn)
      xs = pc.crossSectionNonOriented(E)
      if absorption:
        absorption = NC.createAbsorption(self.fn)
        xs += absorption.crossSectionNonOriented(E)
      if formulaunit:
        xs=xs*sum(self.stoichiometry)
      if plot:
        plt.xscale("log")
        plt.xlabel('Energy [eV]', fontsize="x-large")
        plt.ylabel(' cross section (Barn/formulaunit)', fontsize="x-large")
        plt.grid()
        plt.plot(E, xs, label=self.name)
        plt.legend(fontsize="x-large")
        plt.show()
      if help: print("microscopicCrossSection in Barn/formulaunit as a function of energy in eV")
      return xs
    """MASS ATTENUATION COEFFICIENT (SIGMA/RHO) OF MATERIAL"""
     #_________________________________________________________#
    def massAttenuationCoefficient(self, E, absorption=True, formulaunit=True, plot=False, help=False):
        Na=0.602214076
        massAttenuationCoefficient=self.crossSection(E, absorption, formulaunit)*Na/self.mass
        if plot:
          plt.xscale("log")
          plt.xlabel('Energy [eV]', fontsize="x-large")
          plt.ylabel('Σ/ρ (cm2/g)', fontsize="x-large")
          plt.grid()
          plt.plot(E, massAttenuationCoefficient, label=self.name)
          plt.legend(fontsize="x-large")
          plt.show()
        if help: print("massAttenuationCoefficient in cm2/g as a function of energy in eV")
        return massAttenuationCoefficient
    """MACROSCOPIC CROSS SECTION OF MATERIAL"""
     #__________________________________________#
    def macroscopicCrossSection(self, E, absorption=True, formulaunit=True, plot=False, help=False):
       macroscopicCrossSection=self.massAttenuationCoefficient(E, absorption, formulaunit)*self.density
       if plot:
          plt.xscale("log")
          plt.xlabel('Energy [eV]', fontsize="x-large")
          plt.ylabel('Σ (cm2)', fontsize="x-large")
          plt.grid()
          plt.plot(E, macroscopicCrossSection, label=self.name)
          plt.legend(fontsize="x-large")
          plt.show()
       if help: print("macroscopicCrossSection in cm^2 as a function of energy in eV")
       return macroscopicCrossSection
    """TRANSMISSION OF MATERIAL"""
     #______________________________#
    def transmission(self, E, thickness, absorption=True, formulaunit=True, plot=False, help=False):
        T=np.exp(-self.macroscopicCrossSection(E, absorption, formulaunit)*thickness)*100
        if plot:
          plt.xscale("log")
          plt.xlabel('Energy [eV]', fontsize="x-large")
          plt.ylabel('Transmission (%)', fontsize="x-large")
          plt.grid()
          plt.plot(E, T, label=self.name)
          plt.legend(fontsize="x-large")
          plt.show()
        if help: print("Transmission in % as a function of energy in eV")
        return T
        """THICKNESS OF MATERIAL"""
     #______________________________#
    def findMyThickness(self, Transmission, Energy=1000, absorption=True, formulaunit=True, help=False):
      Z=-np.log(Transmission)/(self.macroscopicCrossSection(Energy, absorption, formulaunit)) #cm
      if help: print("thickness in cm to obtain given transmission")
      return Z
    #def save(self, method, E):
     # s=str(method)
    #  method_name = lambda s: s[s.index('.') + 1 : s.index('(')] if '.' in s and '(' in s else None
     # np.savetxt(f'{self.name}_{method_name(s)}.txt', np.column_stack((E*1000,eval(method))))   # !!!
     # return None




""" Compounds Class"""
class Compound(Material):
    def __init__(self, name, list_of_materials, list_of_weights, density, temperature=293.6, stoichiometry_type="n_of_molecules"):
        print("This class is not finished!!! do not use")
        self.name = name
        formula = ""
        for mat, weight in zip(list_of_materials, list_of_weights):
            formula += mat.formula + str(weight)
        self.formula = formula
        self.density = density
        self.temperature = temperature
        self.list_of_materials = list_of_materials
        self.stoichiometry = stoichiometry_type
        self.list_of_weights = list_of_weights
        SigmaFree=0.
        for mat, weight in zip(self.list_of_materials, self.list_of_weights):
          SigmaFree+=mat.sigmaFree*weight
        self.sigmaFree=SigmaFree
        mass=0.
        for mat, weight in zip(self.list_of_materials, self.list_of_weights):
            mass+=mat.mass*weight
        self.mass=mass

        if stoichiometry_type == "n_of_molecules":
            print("WARNING: list of weights defines the number of molecules, mass percentage can be defined instead by setting stoichiometry_type to 'mass_percentage'")

    """END OF CONSTRUCTOR"""

    def __str__(self):
        return f"\n_____________________\nCOMPOUND PROPERTIES:\n_____________________\n list_of_materials={self.list_of_materials} \n list_of_weights={self.list_of_weights} \n stoichiometry_type={self.stoichiometry_type}\n"

    def crossSection(self, E, absorption=True, formulaunit=True, plot=False, help=False):
      xs=0.
      for mat, weight in zip(self.list_of_materials, self.list_of_weights):
        xs+= mat.crossSection(E, absorption, formulaunit, plot, help)*weight
      if plot:
        plt.xscale("log")
        plt.xlabel('Energy [eV]', fontsize="x-large")
        plt.ylabel(' cross section (Barn/formulaunit)', fontsize="x-large")
        plt.grid()
        plt.plot(E, xs, label=self.name)
        plt.legend(fontsize="x-large")
        plt.show()
      if help: print("microscopicCrossSection in Barn/formulaunit as a function of energy in eV")
      return xs

