# WeakAcidSpeciation
This is a python script for plotting curves of weak acid speciation

You could customize:
1) pkas, up to pka6
2) the total concentration of the weak acid  (default 1)
3) the resolution of the curve  (default 1)
 

10/9/2022 0.3 Updates:
1. Changed WeakAcid into a class with instances which are the dissociation parameters of the object
2. Add a pH parameter in the constructor
3. Print the dissociation parameters of the weak acid at the given pH
4. Remove the costumized step size command line parameter; step size = 0.05 by default
