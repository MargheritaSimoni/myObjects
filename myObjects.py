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
    def __init__(self, name, formula, density, ncmat="AFGA", cif_code=None, functionalGroups=None, temperature=293.6, elements_dictionary=None):
        self.name = name
        self.formula = formula
        self.density = density
        self.temperature = temperature
        self.ncmat = ncmat
        self.cif_code = cif_code
        self.functionalGroups = functionalGroups
        self.elements_dictionary=elements_dictionary
        if not elements_dictionary:
          raise ValueError(f"Elements_dictionary needs to be given as an argument")
        elements, stoichiometry = self.define_Stoichiometry(formula)
        self.elements= elements
        self.stoichiometry =stoichiometry
        for el in self.elements:
          if el not in self.elements_dictionary: raise ValueError(f"Element {el} not found in elements_dictionary")

        self.generate_ncmat()

        """SIGMAFREE OF MATERIAL"""
        SigmaFree=0.
        for el, stoy in zip(self.elements, self.stoichiometry):
          SigmaFree+=elements_dictionary[el].sigmaFree()*stoy
        self.sigmaFree=SigmaFree
        """SIGMABOUND OF MATERIAL"""
        SigmaBound=0.
        for el, stoy in zip(self.elements, self.stoichiometry):
          SigmaBound+=elements_dictionary[el].sigmaBound*stoy
        self.sigmaBound=SigmaBound
        """MASS OF MATERIAL"""
        mass=0.
        for el, stoy in zip(self.elements, self.stoichiometry):
            mass+=elements_dictionary[el].mass*stoy
        self.mass=mass

    """CONSTRUCTOR METHODS"""
            
    """Define material stoichiometry from formula"""
    #______________________________________________#
    def define_Stoichiometry(self,formula):
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

    """build ncmat file"""
    def generate_ncmat(self):
          if self.ncmat=="AFGA":
            if self.functionalGroups==None:
              print("ERROR: When using AFGA model, functional groups must be declared!")
            import subprocess
            command = ['ncrystal_hfg2ncmat','--formula', self.formula,'--spec', self.functionalGroups,'--density', str(self.density),'--title', self.name,'-o', f'{self.name}.ncmat',"--force"]
            creating_ncmat = subprocess.run(command, capture_output=True, text=True)
            self.creating_ncmat= creating_ncmat
            T=str(self.temperature)+"K"
            fn = f'{self.name}.ncmat;temp='+T
          elif self.ncmat=="FreeGas":
            T=str(self.temperature)+"K"
            #fn=f'freegas::'+formula+'/'+str(density)+'gcm3/;temp='+T ######################################################################################### mat1 = NC.createScatter('free_O2.ncmat')
            for el in self.elements:
              c=NC.NCMATComposer(f'freegas::{el}/0.001gcm3/;temp={T}')
              a = c.write(f'free_{el}.ncmat')
              fn="There is an fn for each element, sum is handled inside cross section method"
          elif self.ncmat=="CIF":
            if self.cif_code==None:
              raise ValueError(f"CIF code needs to be given as an argument, as a string: cif_code='codid::1000467'")
            cif_mat = NC.NCMATComposer(self.cif_code)
            a = cif_mat.write(f'{self.name}.ncmat')
            T=str(self.temperature)+"K"
            fn = f'{self.name}.ncmat;temp='+T
          else:
            #print("WARNING: Density and temperature may be defined differently inside NCrystal .ncmat file")
            addT=f";temp=+{self.temperature}"
            fn=self.ncmat+addT
          self.fn=fn

    """END OF CONSTRUCTOR"""
    #______________________________#
    def __str__(self):
        if self.ncmat=="AFGA":
          return f"\n_____________________\nMATERIAL PROPERTIES:\n_____________________\n name={self.name} \n formula={self.formula} \n density={self.density} \n functionalGroups={self.functionalGroups} \n temperature={self.temperature}\n\n This material was created using the Average functional group approximation (AFGA):\n {self.creating_ncmat.stdout}\n {self.creating_ncmat.stderr}\n"
        elif self.ncmat=="FreeGas":
          return f"\n_____________________\nMATERIAL PROPERTIES:\n_____________________\n name={self.name}\n formula={self.formula}\n density={self.density}\n ncmat={self.ncmat}\n temperature={self.temperature}\n\n This material was created using NCrystal FreeGas\n"
        elif self.ncmat=="CIF":
          return f"\n_____________________\nMATERIAL PROPERTIES:\n_____________________\n name={self.name}\n formula={self.formula}\n density={self.density}\n ncmat={self.ncmat}\n temperature={self.temperature}\n\n This material was created from CIF file with code {self.cif_code}\n"
      
        else:
          return f"\n_____________________\nMATERIAL PROPERTIES:\n_____________________\n name={self.name}\n formula={self.formula}\n density={self.density}\n ncmat={self.ncmat}\n temperature={self.temperature}\n\n This material was created using NCrystal libraries\n WARNING: Density and temperature may be defined differently inside NCrystal .ncmat file\n"


    """MICROSCOPIC CROSS SECTION OF MATERIAL"""
     #________________________________________#
    def crossSection(self, E, absorption=True, formulaunit=True, plot=False, help=False):
      if self.ncmat != "FreeGas":
        pc = NC.createScatter(self.fn)
        xs = pc.crossSectionIsotropic(E)
        if absorption:
          absorption = NC.createAbsorption(self.fn)
          xs += absorption.crossSectionIsotropic(E)
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
      if self.ncmat == "FreeGas":
    xs = 0.
    for el, stoy in zip(self.elements, self.stoichiometry):
        scatter_mat = NC.createScatter(f'free_{el}.ncmat')
        absorb_mat = NC.createAbsorption(f'free_{el}.ncmat')  # Assuming this exists
        
        scatter_xs = scatter_mat.crossSectionIsotropic(E)
        absorb_xs = absorb_mat.crossSectionIsotropic(E)
        
        total_xs = scatter_xs + absorb_xs
        xs += total_xs * stoy
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
          plt.ylabel('Σ (cm^-1)', fontsize="x-large")
          plt.grid()
          plt.plot(E, macroscopicCrossSection, label=self.name)
          plt.legend(fontsize="x-large")
          plt.show()
       if help: print("macroscopicCrossSection in cm^2 as a function of energy in eV")
       return macroscopicCrossSection
    """TRANSMISSION OF MATERIAL"""
     #______________________________#
    def transmission(self, E, thickness, percentage=True, absorption=True, formulaunit=True, plot=False, help=False):
        T=np.exp(-self.macroscopicCrossSection(E, absorption, formulaunit)*thickness)
        if percentage:
          T*=100
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
class Compound():
    def __init__(self, name, compound_density, dictionary_of_materials,  dictionary_of_weights, weight_type=None, temperature=293.6):
        self.name = name
        self.compound_density = compound_density
        self.temperature = temperature
        self.dictionary_of_materials = dictionary_of_materials
        self.dictionary_of_weights = dictionary_of_weights
        self.weight_type = weight_type
        self.stoichiometric_weights={}
        self.compute_stoichiometric_weights() #method
        total = sum(self.stoichiometric_weights.values())
        for key in self.stoichiometric_weights:
            self.stoichiometric_weights[key] /= total
        self.sigmaFree=0.
        for name, mat in self.dictionary_of_materials.items():
          print(mat)
          self.sigmaFree+=mat.sigmaFree*self.stoichiometric_weights[name]
        self.mass=0.
        for name, mat in self.dictionary_of_materials.items():
          self.mass+=mat.mass*self.stoichiometric_weights[name]

     # call method inside constructor
    def compute_stoichiometric_weights(self):
        if self.weight_type == "molar":
          self.stoichiometric_weights=self.dictionary_of_weights
        elif self.weight_type == "mass_percentage":
          if sum(self.dictionary_of_weights.values())!=100:
            if sum(self.dictionary_of_weights.values())==1:
              for key in self.dictionary_of_weights:
                self.dictionary_of_weights[key] *= 100
            else: raise ValueError(f"Sum of weights must be 100, not {sum(self.dictionary_of_weights.values())}")
          for mat_name,mat in self.dictionary_of_materials.items():
            self.stoichiometric_weights[mat_name]=self.dictionary_of_weights[mat_name]/mat.mass
        elif self.weight_type == "mass_grams":
          total_mass_g = sum(self.dictionary_of_weights.values())
          for mat_name,mat in self.dictionary_of_materials.items():
            self.stoichiometric_weights[mat_name]=self.dictionary_of_weights[mat_name]/total_mass_g*100/mat.mass
        else:
          raise ValueError(f"Weight type {self.weight_type} not recognized. \n Possible weight types include: \n weight_type='molar' \n weight_type='mass_percentage' \n weight_type='mass_grams'")

    """END OF CONSTRUCTOR"""

    def __str__(self):
        return f"\n_____________________\nCOMPOUND PROPERTIES:\n_____________________\n list_of_materials={self.dictionary_of_materials.keys()} \n list_of_weights={self.dictionary_of_weights} \n weight_type={self.weight_type}\n"

    def crossSection(self, E, absorption=True, formulaunit=True, plot=False, help=False):
      xs=0.
      for name, mat in self.dictionary_of_materials.items():
        xs+= mat.crossSection(E, absorption, True, False, False)*self.stoichiometric_weights[name]

      if formulaunit==False:
        xs=xs/self.sigmaFree
      if plot:
        plt.xscale("log")
        plt.xlabel('Energy [eV]', fontsize="x-large")
        plt.ylabel(' cross section (Barn/formulaunit)', fontsize="x-large")
        plt.grid()
        plt.plot(E, xs, label=self.name)
        plt.legend(fontsize="x-large")
        plt.show()
      if help: print("Stoichiometric ratios of compound are normalized to 1. microscopicCrossSection in Barn/formulaunit as a function of energy in eV")
      return xs

      """MASS ATTENUATION COEFFICIENT (SIGMA/RHO) OF COMPOUND"""
     #_________________________________________________________#
    def massAttenuationCoefficient(self, E, absorption=True, formulaunit=True, plot=False, help=False):
        Na=0.602214076 
        massAttenuationCoefficient=self.crossSection(E, absorption, formulaunit)*Na/(self.mass)
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

    """MACROSCOPIC CROSS SECTION OF COMPOUND"""

    def macroscopicCrossSection(self, E, absorption=True, formulaunit=True, plot=False, help=False):
       Na=0.602214076
       macroscopicCrossSection=self.crossSection(E, absorption, formulaunit)*Na/(self.mass)*self.compound_density
       if plot:
          plt.xscale("log")
          plt.xlabel('Energy [eV]', fontsize="x-large")
          plt.ylabel('Σ (cm^-1)', fontsize="x-large")
          plt.grid()
          plt.plot(E, macroscopicCrossSection, label=self.name)
          plt.legend(fontsize="x-large")
          plt.show()
       if help: print("macroscopicCrossSection in cm^2 as a function of energy in eV")
       return macroscopicCrossSection

    """TRANSMISSION OF COMPOUND"""
     #______________________________#
    def transmission(self, E, thickness, percentage=True, absorption=True, formulaunit=True, plot=False, help=False):
        Na=0.602214076
        T=np.exp(-self.crossSection(E, absorption, formulaunit)*Na/(self.mass)*self.compound_density*thickness)
        if percentage:
          T*=100
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

        """THICKNESS OF COMPOUND"""
        #______________________________#
    def findMyThickness(self, Transmission, Energy=1000, absorption=True, formulaunit=True, help=False):
      Z=-np.log(Transmission)/(self.macroscopicCrossSection(Energy, absorption, formulaunit)) #cm
      if help: print("thickness in cm to obtain given transmission")
      return Z
