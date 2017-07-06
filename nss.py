# -*- coding: utf-8 -*-
"""
Created on Mon Jul 03 09:32:36 2017

@author: gv4723
"""
from __future__ import division
import numpy as np

# import pandas as pd

# NS simple

#def NS(betaV, mats):
#    gam = mats / betaV
#    y = betaV[1] + betaV[2] * ((1 - np.exp(-gam))/(gam)) + betaV[3] * (((1 - np.exp(-gam)) / (gam)) - np.exp(-gam))
#    return  y


# Equation du model de NSS version plus rapide
# nelson--siegel--svensson 2 (a bit faster)
# IN : solution betaV : mats : maturités
# OUT : les taux ou les données....

# Un warning est levé a ce niveau du fait que en python 2 la division par defaut du signe / est entiere...
def NSS2(betaV, mats):
    gam1 = mats / betaV[4]
    gam2 = mats / betaV[5]
    aux1  = 1 - np.exp(-gam1)
    aux2 = 1 - np.exp(-gam2)
    y = betaV[0] + betaV[1] * (aux1 / gam1) + betaV[2] * (aux1 / gam1 + aux1 - 1) + betaV[3] * (aux2 / gam2 + aux2 - 1)
    return (y)

# Fonction objective, la fonction à minimiser qui correspond à la distance entre les données et la fonction
# generic objective function
# In : Vecteur solution optimale betaV, data contient le modele à utiliser pour le calibrage, et les données sur lesquelles calibrer le modèle.
# Out : distance entre modèle et données réelles

def OF(betaV, data):

    mats =data.get("mats") # Les maturités
    yM = data.get("yM") # Les taux
    # model = data.get("model") # Le nom du modèle à utiliser, NS ou NSS2
    y = NSS2(betaV,mats)

    # crossprod() : Somme des carrées de ,
    # qui correspond donc à la distance entre le modèle y et les données yM
    # Pour changer la distance utilisé, c'est ici que ca se passe :
    aux=(y-yM)
    aux =  np.dot(aux,aux)
    return(aux)

#################### FONCTIONS CONCUES POUR ETRE UTILISE AVEC DES ARRAYS

# Fontions concues pour etre utilisée en DE
def mRU(m,n):
    return np.random.uniform(size=(m,n))

def mRN(m,n):
    return np.random.normal(size=(m,n))

def shift(x):
    rr = len(x)
    return np.hstack((x[rr-1],np.array(x[range(rr-1)])))

def soustrac1(x, maxV):
    return x - maxV

def pen(mP, pso, vF):
    minV = pso.get("min")
    maxV = pso.get("max")
    ww = pso.get("ww")
    
    # max constraint: if larger than maxV, element in A is positiv
    # A = mP - (maxV)
    A = np.apply_along_axis(soustrac1,0,mP,maxV=maxV)
    A = A + abs(A)

    # max constraint: if smaller than minV, element in B is positiv
    # B = minV - mP
    B = np.apply_along_axis(soustrac1,0,-mP,maxV=-maxV)
    B = B + abs(B)

    ## beta 1 + beta2 > 0
    C = ww * ((mP[0,:] + mP[1,:]) - abs(mP[0,:] + mP[1,:]))
    A = ww * np.sum((A + B),axis = 0) * vF - C
    return(A)


def DE(de,dataList,OF):

    # main algorithm
    # ------------------------------------------------------------------
    # set up initial population
    #mP = de.get("min") + np.dot(np.diag(de.get("max") - de.get("min")) , mRU(de.get("d"),de.get("nP")))
    mP = np.apply_along_axis(soustrac1, 0, np.dot(np.diag(de.get("max") - de.get("min")) , mRU(de.get("d"),de.get("nP"))) ,maxV = -de.get("min"))

    # include extremes
    mP[:,0:de.get("d")] = np.diag(de.get("max"))
    mP[:,(de.get("d")):(2*de.get("d"))] = np.diag(de.get("min"))

    # evaluate initial population
    vF = np.apply_along_axis(OF, 0, mP,data = dataList)
    #vF = apply(mP,2,OF,data = dataList)

    # constraints
    vP = pen(mP,de,vF)
    vF = vF + vP

    # keep track of OF
    Fmat = np.empty((de.get("nG"), de.get("nP")), dtype=float)

    for g in range(de.get("nG")):
        # update population
        vI = np.random.choice(range(de.get("nP")), de.get("nP"))
        R1 = shift(vI)
        R2 = shift(R1)
        R3 = shift(R2)
        
        # prelim. update
        mPv = mP[:,R1] + de.get("F") * (mP[:,R2] - mP[:,R3])
        
        if de.get("R") > 0:
            mPv = mPv + de.get("R") * mRN(de.get("d"),de.get("nP"))
        
        #
        mI = mRU(de.get("d"),de.get("nP")) > de.get("CR")
        mPv[mI] = mP[mI]
        
        #
        vFv = np.apply_along_axis(OF, 0, mPv,data = dataList)
        
        # constraints
        vPv = pen(mPv,de,vF)
        vFv = vFv + vPv
        vFv[(np.isfinite(vFv))] = 1000000
        
        # find improvements
        logik = vFv < vF
        mP[logik] = mPv[logik]
        vF[logik] = vFv[logik]
        Fmat = vF

    # g in 1:nG
    sGbest = np.min(vF)
    sgbest = np.argmin(vF)

    # return best solution
    return ({"beta" : mP[:,sgbest], "OFvalue" : sGbest, "popF" : vF, "Fmat" : Fmat})


# application de l'evolution differentielle appliqué à Nelson siegel svennsson
# maturites 
mats = np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.])

# taux
yM = np.array([0.1,0.2,0.3,0.4,0.5,0.55,0.57,0.59,0.61])

# 
de = {"min"	: np.array([0.,-15.,-30.,-30.,0.  ,2.5]),# minimum du vecteur des paramètres à estimer 
            "max"	: np.array([15., 30., 30., 30.,2.5,5. ]),# Le maximum du vecteur des paramètres
		"d"	: 6, # Dimension du param-tre à estimer ou effectif de la population
		"nP"	: 200,
		"nG"	: 600,
		"ww"	: 0.1,
		"F"	: 0.50,
		"CR"	: 0.99,
		"R"	: 0.00 # random term (added to change)
}

dataList = {"yM" : yM, "mats" : mats, "model" : NSS2}

sols = DE(de = de,dataList = dataList,OF = OF)
