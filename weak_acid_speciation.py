'''
Weak Acid Speciation Plotter
weak_acid_speciation.py
0.3 edition
Modified: 10/9/2022

10/9/2022 0.3 Updates:
1. Changed WeakAcid into a class with instances which are the dissociation parameters of the object
2. Add a pH parameter in the constructor
3. Print the dissociation parameters of the weak acid at the given pH
4. Remove the costumized step size command line parameter; step size = 0.05 by default
'''




import math as m
import numpy as np
import matplotlib.pyplot as plt
import sys



# a class that holds an weak acid instance with parameters
class WeakAcid:

    #contructor that intializes instance variables
    def __init__(self, pKa_list, pH, TA = 1):
        self.pH = pH
        self.H_Conc = self.pH_to_H(pH)
        self.pKaList = pKa_list
        self.KaList = self.pka_to_Ka(pKa_list)
        self.TA = TA

        #buffer capacity
        self.bufferCapacity = self.buffer_capacity(self.KaList, self.H_Conc, self.TA)

        #degree of dissociation
        self.degreeOfDissociation = self.degree_of_dissociation(self.KaList, self.H_Conc)
        
        #alpha values
        self.alpha0 = self.alpha(self.KaList, self.H_Conc, 0)
        self.alpha1 = self.alpha(self.KaList, self.H_Conc, 1)
        self.alpha2 = self.alpha(self.KaList, self.H_Conc, 2)
        self.alpha3 = self.alpha(self.KaList, self.H_Conc, 3)
        self.alpha4 = self.alpha(self.KaList, self.H_Conc, 4)
        self.alpha5 = self.alpha(self.KaList, self.H_Conc, 5)

        #resolution of the graph
        self.stepSize = 0.05
        

    """Convert pKa_list"""
    def pka_to_Ka(self, pKa_list):
        # if the list has greater than 6 elements or empty, terminate and print
        if len(pKa_list) == 0 or len(pKa_list) > 6:
            print("pKa values invalid")
            exit()

        #there are 6 pKa in total; program automatically identifies the length of pka list and then append 35 to the rest 6-len(pKa_list) elememts
        for i in range(6-len(pKa_list)):
            pKa_list.append(100)

        #create Ka_list according to pKa_list
        Ka_list = []
        for i in pKa_list:
            Ka_list.append(10**(-i))
        
        #return the Ka list
        return Ka_list



    '''return the H+ concentration at given pH'''
    def pH_to_H (self, pH):
        H_conc = 10**(-pH)
        return H_conc



    '''return the list H+ concentration given pH list'''
    def pH_to_H_list (self, pH_list):
        H_conc_list = []
        for i in pH_list:
            H_conc_list.append(self.pH_to_H(i))
        return H_conc_list



    """returns the denominator given a Ka_list and the H+ concentration"""
    def denominator(self, H_conc):
        denom = H_conc**6 + (H_conc**5)*self.KaList[0] + (H_conc**4)*self.KaList[0]*self.KaList[1] + (H_conc**3)*self.KaList[0]*self.KaList[1]*self.KaList[2] + (H_conc**2)*self.KaList[0]*self.KaList[1]*self.KaList[2]*self.KaList[3] + (H_conc)*self.KaList[0]*self.KaList[1]*self.KaList[2]*self.KaList[3]*self.KaList[4] + self.KaList[0]*self.KaList[1]*self.KaList[2]*self.KaList[3]*self.KaList[4]*self.KaList[5]
        return denom



    """returns the alpha values based on Ka list and H+ concentration, which indicates which alpha"""
    def alpha(self, Ka_list, H_conc, which):
        # call the denominator function to calculate the denominator
        denom = self.denominator(H_conc)

        # calculate alpha concentration based on that dummy equation
        alpha = (H_conc**(6-which)) * np.prod(Ka_list[:which]) / denom
        return alpha



    '''return the degree of dissociation by given Ka_list and H_con'''
    def degree_of_dissociation(self, Ka_list, H_conc):
        d = 0
        # calculate degree of dissociation based on equation
        for i in range(6):
            d +=  i * self.alpha(Ka_list, H_conc, i)
        return d



    '''return the buffer capacity based on given total weak acid concentration and the alpha values'''
    def buffer_capacity(self, Ka_list, H_conc, Total_Acid_Conc):
        #calculated the buffer capacity accoridng to equation (the sum of neighbour products times the total weak acid concentration)
        sum = 0
        for k in range(5):
            sum += self.alpha(Ka_list, H_conc, k) * self.alpha(Ka_list, H_conc, k+1)

        return sum * Total_Acid_Conc



    "a function that prints all the parameters"
    def printAll(self):
        print("Given: pH =", self.pH, "  Given pKas:", self.pKaList, "  TA:", self.TA, "Moles")
        print("α0 = ", self.alpha0)
        print("α1 = ", self.alpha1)
        print("α2 = ", self.alpha2)
        print("α3 = ", self.alpha3)
        print("α4 = ", self.alpha4)
        print("α5 = ", self.alpha5)
        print("Buffer Capacity = ", self.bufferCapacity)
        print("Degree of Dissociation = ", self.degreeOfDissociation)



    '''plot the alphas, degree of dissociation, and buffer capacity '''
    '''return alpha lists, degree of dissociation list, and buffer capacity list'''
    def plotter(self):
        #create list of pH
        pH_list = np.arange(0.0, 14.0 + self.stepSize, self.stepSize)

        #create list of H_conc
        H_conc_list = self.pH_to_H_list(pH_list)

        #create empty list of alpha values list
        alphas_list = []

        #create lists for alphas
        for i in range(6):
            tempList = []
            for H_Conc in H_conc_list:
                tempList.append(self.alpha(self.KaList, H_Conc, i))
            alphas_list.append(tempList)

        #create list for weak acid dissociation and buffer) capacity
        degree_dissociation_list = []
        buffer_capacity_list = []
        for H_Conc in H_conc_list:
            degree_dissociation_list.append(self.degree_of_dissociation(self.KaList, H_Conc))
            buffer_capacity_list.append(self.buffer_capacity(self.KaList, H_Conc, self.TA))

        #plot the degree of dissociation and buffere capacity
        plt.plot(pH_list, degree_dissociation_list, label = 'ŋ')
        plt.plot(pH_list, buffer_capacity_list, label = 'β')

        #create a list for labels of lines
        line_label_list = ["α0", "α1", "α2", "α3", "α4", "α5", "α6"]

        #use for loop to plot the alphas with pH, labeling them with the string list just created
        for i in range(len(list(sys.argv[1].split(',')))+1):
            plt.plot(pH_list, alphas_list[i], label = line_label_list[i])  #alphas_list contains all the sublists of alphas

        #insert a line showing the specified pH
        if len(sys.argv) > 2:
            plt.axvline(x = float(sys.argv[2]), color = 'black', label = "pH=" + str(sys.argv[2]))
        
        #create a title
        plt.title(label = "Speciation Curve of Acid with pKas: " + " [" + str(sys.argv[1]) + "]")

        #show the legends and plot
        plt.legend()
        plt.show()

        #return the lists created to be reused in other cases
        return alphas_list, degree_dissociation_list, buffer_capacity_list

        


    '''print user statements'''
    def printUsage():
        print('')
        print('Instructions:')
        print('1. Check if you set the working directory to the correct source file directory')
        print('2. Only need to enter eligible pKas, the rest of the pKas will be auto-filled', '\n')
        print('Usages')
        print('Usage(Windows): py weak_acid_speciation.py <pKas separated by commas> <pH> <total acid concentration in Moles = 1>')
        print('E.G. (Windows): py weak_acid_speciation.py 2.148,7.198,12.375 7.198 1')  #the example is phosphate
        print('Usage(Mac): python3 weak_acid_speciation.py <pKas separated by commas> <pH> <total acid concentration in Moles = 1>')
        print('E.G. (Mac): python3 weak_acid_speciation.py 2.148,7.198,12.375 7.198 1')
        exit()



#create the main function that uses command line arguments to control plots
if __name__ == '__main__':

    if len(sys.argv) > 3:
        acid = WeakAcid([float(i) for i in list(sys.argv[1].split(','))], float(sys.argv[2]), float(sys.argv[3]))
        
    if len(sys.argv) == 3: 
        acid = WeakAcid([float(i) for i in list(sys.argv[1].split(','))], float(sys.argv[2]), TA = 1)
        

    if len(sys.argv) == 2: 
        acid = WeakAcid([float(i) for i in list(sys.argv[1].split(','))], pH = 7, TA = 1)


    #print the acid dissociation parameters at given pH
    acid.printAll()

    #plot the dissociation curve
    acid.plotter()


    #print out the usage instruction of the program if users enter invalid arguments
    if len(sys.argv) < 2:
        acid.printUsage()