#Weak acid speciation plotter

import math as m
import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import AutoMinorLocator



class WeakAcid:
    pH = 0
    H_conc = 0
    Ka_list = []
    Total_Acid_Concentration = 1

    def __init__(self, pKa_list, pH, Total_Acid_Concentration):
        WeakAcid.pH = pH
        WeakAcid.H_conc = self.pH_to_H(pH)
        WeakAcid.Ka_list = self.pka_to_Ka(pKa_list)
        WeakAcid.Total_Acid_Concentration = Total_Acid_Concentration


    """Convert pKa_list"""
    def pka_to_Ka(pKa_list):
        # if the list has greater than 6 elements or empty, terminate and print
        if len(pKa_list) == 0 or len(pKa_list) > 6:
            print("pKa values invalid")
            exit()

        #there are 6 pKa in total; program automatically identifies the length of pka list and then append 35 to the rest 6-len(pKa_list) elememts
        for i in range(6-len(pKa_list)):
            pKa_list.append(50)

        #create Ka_list according to pKa_list
        Ka_list = []
        for i in pKa_list:
            Ka_list.append(10**(-i))
        
        #return the Ka list
        return Ka_list



    '''return the H+ concentration at given pH'''
    def pH_to_H (pH):
        H_conc = 10**(-pH)
        return H_conc



    '''return the list H+ concentration given pH list'''
    def pH_to_H_list (pH_list):
        H_conc_list = []
        for i in pH_list:
            H_conc_list.append(WeakAcid.pH_to_H(i))
        return H_conc_list



    """returns the denominator given a Ka_list and the H+ concentration"""
    def denominator(Ka_list, H_conc):
        denom = H_conc**6 + (H_conc**5)*Ka_list[0] + (H_conc**4)*Ka_list[0]*Ka_list[1] + (H_conc**3)*Ka_list[0]*Ka_list[1]*Ka_list[2] + (H_conc**2)*Ka_list[0]*Ka_list[1]*Ka_list[2]*Ka_list[3] + (H_conc)*Ka_list[0]*Ka_list[1]*Ka_list[2]*Ka_list[3]*Ka_list[4] + Ka_list[0]*Ka_list[1]*Ka_list[2]*Ka_list[3]*Ka_list[4]*Ka_list[5]
        return denom



    """returns the alpha values based on Ka list and H+ concentration, which indicates which alpha"""
    def alpha(Ka_list, H_conc, which):
        # call the denominator function to calculate the denominator
        denom = WeakAcid.denominator(Ka_list, H_conc)

        # calculate alpha concentration based on that dummy equation
        alpha = (H_conc**(6-which)) * np.prod(Ka_list[:which]) / denom
        return alpha



    '''return the degree of dissociation by given Ka_list and H_con'''
    def degree_of_dissociation(Ka_list, H_conc):
        d = 0
        # calculate degree of dissociation based on equation
        for i in range(6):
            d +=  i * WeakAcid.alpha(Ka_list, H_conc, i)
        return d



    '''return the buffer capacity based on given total weak acid concentration and the alpha values'''
    def buffer_capacity(Ka_list, H_conc, Total_Acid_Conc):
        #calculated the buffer capacity accoridng to equation (the sum of neighbour products times the total weak acid concentration)
        sum = 0
        for k in range(5):
            sum += WeakAcid.alpha(Ka_list, H_conc, k) * WeakAcid.alpha(Ka_list, H_conc, k+1)

        return sum * Total_Acid_Conc



    '''plot the alphas, degree of dissociation, and buffer capacity '''
    def plotter(pKa_list, Total_Acid_Conc, step_size):
        #get a list of Ka
        Ka_list = WeakAcid.pka_to_Ka(pKa_list)

        #create list of pH
        pH_list = np.arange(0.0, 14.0 + step_size, step_size)

        #create list of H_conc
        H_conc_list = WeakAcid.pH_to_H_list(pH_list)

        #create empty list of alpha values list
        alphas_list = []

        #create lists for alphas
        for i in range(6):
            tempList = []
            for H_Conc in H_conc_list:
                tempList.append(WeakAcid.alpha(Ka_list, H_Conc, i))
            alphas_list.append(tempList)

        #create list for weak acid dissociation and buffer) capacity
        degree_dissociation_list = []
        buffer_capacity_list = []
        for H_Conc in H_conc_list:
            degree_dissociation_list.append(WeakAcid.degree_of_dissociation(Ka_list, H_Conc))
            buffer_capacity_list.append(WeakAcid.buffer_capacity(Ka_list, H_Conc, Total_Acid_Conc))

        #plot the degree of dissociation and buffere capacity
        plt.plot(pH_list, degree_dissociation_list, label = 'ŋ')
        plt.plot(pH_list, buffer_capacity_list, label = 'β')

        #create a list for labels of lines
        line_label_list = ["α0", "α1", "α2", "α3", "α4", "α5", "α6"]

        #use for loop to plot the alphas with pH, labeling them with the string list just created
        for i in range(len(list(sys.argv[1].split(',')))+1):
            plt.plot(pH_list, alphas_list[i], label = line_label_list[i])  #alphas_list contains all the sublists of alphas


        #show the legends
        plt.legend()
        plt.show()
        


#create the main function that use command line arguments that allow control by terminal
if __name__ == '__main__':
    if len(sys.argv)>3:
        pKalist = [float(i) for i in list(sys.argv[1].split(','))]
        Total_Acid_Conc = float(sys.argv[2])
        step_size = float(sys.argv[3])
        WeakAcid.plotter(pKalist, Total_Acid_Conc, step_size)

    elif len(sys.argv) == 3:
        pKalist = [float(i) for i in list(sys.argv[1].split(','))]
        Total_Acid_Conc = float(sys.argv[2])
        step_size = 0.05
        WeakAcid.plotter(pKalist, Total_Acid_Conc, step_size)


    elif len(sys.argv) == 2:
        pKalist = [float(i) for i in list(sys.argv[1].split(','))]
        Total_Acid_Conc = 1
        step_size = 0.05
        WeakAcid.plotter(pKalist, Total_Acid_Conc, step_size)
    

    #print out the usage instruction of the program if users enter invalid arguments
    else:
        print('')
        print('Instructions:')
        print('1. Check if you set the working directory to the correct source file directory')
        print('2. Only need to enter eligible pKas, the rest of the pKas will be auto-filled', '\n')
        print('Usages')
        print('Usage(Windows): py weak_acid_speciation.py <pKas separated by commas> <total acid concentration in Moles = 1> <stepsize = 0.05>')
        print('E.G. (Windows): py weak_acid_speciation.py 2.148,7.198,12.375 1 0.05')  #the example is phosphate
        print('Usage(Mac): python3 weak_acid_speciation.py <pKas separated by commas> <total acid concentration in Moles = 1> <stepsize = 0.05>')
        print('E.G. (Mac): python3 weak_acid_speciation.py 2.148,7.198,12.375 1 0.05')
        exit()
