# -*- coding: utf-8 -*-

"""
Created on Sun Mar  5 13:41:37 2017

@author: workenv
"""

# CALIBRATION DU MODEL DE NILSON  SIEGEL SVENSSON

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc


df = pd.read_csv("yield_curve_spot_bce.csv", header=None, decimal=",")

to=df[1]

# Fonction d'optimisatiion, 
# minimisation de la somme des écarts entre taux modélisés et taux théoriques 

def fun(x,to):
    return sum(((x[0]+x[1]*((1-np.exp(-to/x[4]))/(to/x[4]))+x[2]*(((1-np.exp(-to/x[4]))/(to/x[4]))-np.exp(-to/x[4]))+x[3]*(((1-np.exp(-to/x[5]))/(to/x[5]))-np.exp(-to/x[5])))-df[1])**2)

# Optimisation sans contraintes, méthode du simplexe

res_sc =scipy.optimize.minimize(fun, (.1,1,1,1,1,1),args=(to), method= "Nelder-Mead")

# optimisation sous contraintes

cons = ({'type': 'ineq', 'fun': lambda x:  x[1] + x[2]})

bnds = ((0, None), (None, None),(None, None), (None, None),(0, None), (0, None))

res_c =sc.optimize.minimize(fun, (0.1,1,1,1,1,1),args=(to), method= "SLSQP",bounds=bnds)

#Fonction de calcul des courbes modélisées
# ycm = yields curves estimated

def ycm(to, x):
    return x[0]+x[1]*((1-np.exp(-to/x[4]))/(to/x[4]))+x[2]*(((1-np.exp(-to/x[4]))/(to/x[4]))-np.exp(-to/x[4]))+x[3]*(((1-np.exp(-to/x[5]))/(to/x[5]))-np.exp(-to/x[5]))
    
    #Graphiqques 

plt.plot(df[0], df[1] ,'o-', df[0], ycm(to,res_c.x),'g-')
plt.title('BCE Yield Curves')
plt.ylabel('Rates')
plt.ylabel('Maturity')

#f2 = sc.interpolate.interp1d(df[0], df[1], kind='cubic')
#plt.plot(df[0], df[1] ,'o-', df[0], f2(df[0]),'g-')
#
#tck = sc.interpolate.splrep(df[0], df[1], s=0)

sc.interp