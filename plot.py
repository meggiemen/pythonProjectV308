import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import sqrt
import pandas as pd
import scipy.constants as const
from scipy.optimize import curve_fit                        # Funktionsfit:     popt, pcov = curve_fit(func, xdata, ydata) 
from uncertainties import ufloat                            # Fehler:           fehlerwert =  ulfaot(x, err)
from uncertainties import unumpy as unp 
from uncertainties.unumpy import uarray                     # Array von Fehler: fehlerarray =  uarray(array, errarray)
from uncertainties.unumpy import (nominal_values as noms,   # Wert:             noms(fehlerwert) = x
                                  std_devs as stds)         # Abweichung:       stds(fehlerarray) = errarray
from math import isclose
  
# Definiton und Funktion des relativen Prozentualen Fehlers Zwischen X_app und X_theo
def rel_err(X_app,X_theo):
    return (np.abs(X_app-X_theo)/X_theo)*100

### Einzelspulen Plot

###Einzelspule B-Feld
def B_Spule(lPosition,lWindungen,lRadius,lStrom):
    return (const.mu_0*lStrom*lWindungen*lRadius**2)/(2*((lRadius**2+lPosition**2)**(3/2)))

x_Einzel1, B_Einzel1 = np.genfromtxt('tables/EinzelSpule1.txt',unpack=True,skip_header=1) #Messdaten importieren
x_Einzel2, B_Einzel2 = np.genfromtxt('tables/EinzelSpule2.txt',unpack=True,skip_header=1) 

###Erste Spule

plt.plot(x_Einzel1,B_Einzel1,'.')
plt.xlabel(r'$x\,/\,\unit{cm}$')
plt.ylabel(r'$B\,/\,\unit{mT}$')

plt.savefig('build/EinzelspulePlot1.pdf',bbox_inches='tight')
plt.clf()

###Zweite Spule

plt.plot(x_Einzel2,B_Einzel2,'.')
plt.xlabel(r'$x\,/\,\unit{cm}$')
plt.ylabel(r'$B\,/\,\unit{mT}$')

plt.savefig('build/EinzelspulePlot2.pdf',bbox_inches='tight')
plt.clf()

### HelmholtzSpulenPaar

I_Helm=3
N_Helm=100
R_Helm=0.0625
breite_Helm=0.033

###Koordinatenverschiebung
def x_Neu(lAlt,lDistanz):
    return lAlt-lDistanz/2

###Magnetfeld Theorie berechnungsformel
def B_helm(lPosition,lAbstand):
    return B_Spule(lPosition-(lAbstand/2),N_Helm,R_Helm,I_Helm)+B_Spule(lPosition+(lAbstand/2),N_Helm,R_Helm,I_Helm)

###Distanz 1

x_Helmholtz1, B_Helmholtz1 = np.genfromtxt('tables/SpulenPaar7.txt',unpack=True,skip_header=1)
x_Helmholtz1=x_Neu(x_Helmholtz1,7)

x1=np.linspace(x_Helmholtz1[0],x_Helmholtz1[-1],100)
B1=B_helm(x1/100,0.07)*10**(3)

plt.plot(x1,B1,'-k',label='Theorie')
plt.plot(x_Helmholtz1,B_Helmholtz1,'.',label='Messwerte')
plt.xlabel(r'$x\,/\,\unit{cm}$')
plt.ylabel(r'$B\,/\,\unit{mT}$')
plt.legend(loc='best')

plt.savefig('build/HelmholtzPlot1.pdf',bbox_inches='tight')
plt.clf()

###Distanz 2

x_Helmholtz2, B_Helmholtz2 = np.genfromtxt('tables/SpulenPaar12.txt',unpack=True,skip_header=1)
x_Helmholtz2=x_Neu(x_Helmholtz2,12)

x2=np.linspace(x_Helmholtz2[0],x_Helmholtz2[-1],100)
B2=B_helm(x2/100,0.12)*10**(3)

plt.plot(x2,B2,'-k',label='Theorie')
plt.plot(x_Helmholtz2,B_Helmholtz2,'.',label='Messwerte')
plt.xlabel(r'$x\,/\,\unit{cm}$')
plt.ylabel(r'$B\,/\,\unit{mT}$')
plt.legend(loc='best')

plt.savefig('build/HelmholtzPlot2.pdf',bbox_inches='tight')
plt.clf()

###Distanz 3

x_Helmholtz3, B_Helmholtz3 = np.genfromtxt('tables/SpulenPaar18.txt',unpack=True,skip_header=1)
x_Helmholtz3=x_Neu(x_Helmholtz3,18)

x3=np.linspace(x_Helmholtz3[0],x_Helmholtz3[-1],100)
B3=B_helm(x3/100,0.18)*10**(3)

plt.plot(x3,B3,'-k',label='Theorie')
plt.plot(x_Helmholtz3,B_Helmholtz3,'.',label='Messwerte')
plt.xlabel(r'$x\,/\,\unit{cm}$')
plt.ylabel(r'$B\,/\,\unit{mT}$')
plt.legend(loc='best')

plt.savefig('build/HelmholtzPlot3.pdf',bbox_inches='tight')
plt.clf()

### HystereseKurve

###Import
I_HystereseNeu,B_HystereseNeu = np.genfromtxt('tables/HystereseNeukurve.txt',unpack=True,skip_header=1)
I_HystereseMinus,B_HystereseMinus,I_HysteresePlus,B_HysteresePlus = np.genfromtxt('tables/HystereseKurve.txt',unpack=True,skip_header=1)

###Korrektur für die Feldstärke
H_HystereseNeu = I_HystereseNeu*595/(2*np.pi*0.135)
H_HystereseMinus = I_HystereseMinus*595/(2*np.pi*0.135)
H_HysteresePlus = I_HysteresePlus*595/(2*np.pi*0.135)

###Fit an Neukurve und Permeabilität
# def Neu_fit(a,b,c,x):
#     return a*np.exp(b*x)+c

# a_guess = -0.6
# b_guess = -0.02
# c_guess = 0.8   

# popt,pcov = curve_fit(Neu_fit,H_HystereseNeu/1000,B_HystereseNeu/1000,p0=(a_guess,b_guess,c_guess))
# Neu_H_Werte = np.linspace(H_HystereseNeu[0]/1000,H_HystereseNeu[-1]/1000,100)

# print(popt)

# plt.plot(Neu_H_Werte,Neu_fit(popt[0],-popt[1],popt[2],Neu_H_Werte),'k-',label=r'Fit der Form $a\cdot e^{b\cdot x+c}+d$')
# plt.plot(H_HystereseNeu/1000,B_HystereseNeu/1000,'r--',label='Neukurve')

# plt.legend(loc='best')
# plt.grid(visible=True)
# plt.xlabel(r'$H / \frac{\unit{A}}{\unit{m}}$')
# plt.ylabel(r'$B\,/\,\unit{mT}$')

# plt.savefig('build/HystereseFitPlot.pdf',bbox_inches='tight')
# plt.clf()

###Berechnung der Koerzitifkraft und Remanenz
x1=np.linspace(H_HystereseMinus[18],H_HystereseMinus[23],100)
x2=np.linspace(-464.26803,-464.26805,100000)
R_fit = np.polyfit(H_HystereseMinus[18:23],B_HystereseMinus[18:23],1)
K_fit = np.polyfit(H_HystereseMinus[20:23],B_HystereseMinus[20:23],1)

R_fit = np.poly1d(R_fit)
K_fit = np.poly1d(K_fit)

for i in x2:
    # if K_fit(i)<0.000000001 and K_fit(i)>-0.000000001:
    #     print(K_fit(i),i)
    if isclose(K_fit(i),2.4840574042173102e-11): 
        print('Koerzitivkraft:',i)
        

print('Remanenz:',R_fit(0))

##Plot der Neukurve
plt.plot(H_HystereseNeu,B_HystereseNeu,'r--',label='Neukurve')

##Plot der  Hysteresekurve
plt.plot(H_HysteresePlus,B_HysteresePlus,'b--',label='Hysteresekurve')
plt.plot(H_HystereseMinus,B_HystereseMinus,'b--')

##Plot der Sättigung
plt.plot(H_HystereseMinus[0],B_HystereseMinus[0],'.',color='k',label='Sättigung')
plt.plot(H_HysteresePlus[0],B_HysteresePlus[0],'.',color='k')

##Plot der Remanenz
plt.plot(H_HystereseMinus[20],B_HystereseMinus[20],'.',color='g',label='Remanenz')
plt.plot(H_HysteresePlus[20],B_HysteresePlus[20],'.',color='g')

##Plot der Koerzitifkraft
plt.plot(H_HystereseMinus[21]-70,0,'.m',label='Koerzitivkraft')
plt.plot(H_HysteresePlus[21],B_HysteresePlus[21],marker='.',color='m')

plt.legend(loc='best')
plt.grid(visible=True)
plt.xlabel(r'$H / \frac{\unit{A}}{\unit{m}}$')
plt.ylabel(r'$B\,/\,\unit{mT}$')

plt.savefig('build/HysteresePlot.pdf',bbox_inches='tight')
plt.clf()
