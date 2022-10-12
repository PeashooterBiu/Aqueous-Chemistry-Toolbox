weak_acid_speciation.py
This is a toolbox for plotting curves of weak acid speciation given a list of pKas

You could customize:
1) pkas, up to pka6
2) the total concentration of the weak acid  (default 1)
3) the resolution of the curve  (default 1)
 
10/9/2022 Updates:
1. Changed WeakAcid into a class with instances which are the dissociation parameters of the object
2. Add a pH parameter in the constructor
3. Print the dissociation parameters of the weak acid at the given pH
4. Remove the costumized step size command line parameter; step size = 0.05 by default



solubilitySimplified.py
This is a toolbox for calculating solubility interms of cation and anion given pKa lists of the cation and the anion and a pH.

10/12/2022 Updates:
1. Show a plot of solubilities of cation and anion against pH
