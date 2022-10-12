'''
Solublility plotter & calculator for general salts 
Solubility.py
0.1 edition
Modified: 10/9/2022
'''


import weak_acid_speciation as wa
import activity as a
import numpy as np
import matplotlib.pyplot as plt


'''
    holds Cation and Anion, both weak acids
    if one of them is not weak acid, the method still works as long as the pka is correct (very negative)
'''
class SimpleSolubility:

    '''constructor that creates cation and anion as instances of the WeakAcid class'''
    #μ0 is the initial ionic strenght of the system before the salt dissolves
    #nCation and nAnion are the coefficients in the salt
    def __init__(self, pKasCation, pKasAnion, Ksp, pH, nCation, nAnion, whichCation, whichAnion):
        self.Cation = wa.WeakAcid(pKasCation, pH)
        self.Anion = wa.WeakAcid(pKasAnion, pH)
        self.pH = pH
        self.Hconc = self.Cation.pH_to_H(self.pH)
        self.Ksp = Ksp
        self.whichCation = whichCation
        self.whichAnion = whichAnion
        self.nCation = nCation
        self.nAnion = nAnion
        self.CationConc, self.AnionConc = self.eqConc(self.Hconc)
        self.KspPrime = self.Ksp_prime(self.Hconc)


    #calculate and return the Ksp_prime
    def Ksp_prime(self, Hconc):
        Ksp_prime =  self.Ksp / (self.Cation.alpha(self.Cation.KaList, Hconc, self.whichCation))**self.nCation / (self.Anion.alpha(self.Anion.KaList, Hconc, self.whichAnion))**self.nAnion
        return Ksp_prime


    #calculate and return the equilirbium concentration of anion and cation
    def eqConc(self, Hconc):
        Ksp_prime = self.Ksp_prime(Hconc)
        AnionConc = ((self.nAnion / self.nCation) ** self.nCation * Ksp_prime) ** (1 / (self.nAnion + self.nCation))
        CationConc = ((self.nCation / self.nAnion) ** self.nAnion * Ksp_prime) ** (1 / (self.nAnion + self.nCation))
        return CationConc, AnionConc
    

    def printAll(self):
        print("Given pH =", self.pH, ", the solubilities are:")
        print("[Cation] =", self.CationConc, "Mole")
        print("[Anion] =", self.AnionConc, "Mole")
        print("Ksp_prime = ", self.KspPrime)


    #plotter create discrete value lists to be plotted
    def plotter(self):
        #create list of pH
        pH_list = np.arange(0.0, 14.0 + 0.05, 0.05)

        #create list of H_conc
        H_conc_list = self.Cation.pH_to_H_list(pH_list)

        #create lists for solubility of Cation and Anion
        CationList = []
        AnionList = []

        for H_conc in H_conc_list:
            CationConc, AnionConc = self.eqConc(H_conc)
            CationList.append(float(CationConc))
            AnionList.append(float(AnionConc))


        plt.plot(pH_list, CationList, label = "[Cation]")
        plt.plot(pH_list, AnionList, label = "[Anion]")
        plt.ylabel("Solubility(moles)")
        plt.xlabel("pH")

        plt.legend()
        plt.show()



if __name__ == '__main__':
    calciumPhosphate = SimpleSolubility([7.7,9.4,11], [2.148, 7.198, 12.375], 2e-28, 7, 3, 2, 0, 3)
    calciumPhosphate.printAll()

    calciumPhosphate.plotter()