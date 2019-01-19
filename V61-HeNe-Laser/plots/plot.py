import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 1000, 1000)
mhub = const.value('Bohr magneton') # das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt


def stabiparams(L,r1,r2):
	return ((L**2)/(r1*r2) -L*((1/r1) + (1/r2)) + 1)


def bed(L, r1):
     return (1 - L / r1)

def mittel(x):              #the real mean()-ing of life
	return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))


def relf(l,m):  #in Prozent
	return (np.absolute(l-m)/l)*100

######Vorbereitungsaufgabe######
## r Parameter der Spiegel
roc = 140 # in centimeter out coupling 
rk1 = 100 #konkav 
rk2 = roc 
L1 = np.linspace(-10, 290, 10**3)
L2 = np.linspace(-10, 150, 10**3)

plt.plot(L1, stabiparams(L1,roc,rk2), label=r'$r_1 = 1400$mm & $r_2= 1400$mm')
plt.xlabel(r'Resonatorabstand $L/cm$')
plt.ylabel(r'Stabilitätsparameter $g_1 \cdot g_2$')
plt.axhline(1, linestyle = '--',color='r', label='Laserbereich')
plt.axhline(0, linestyle = '--',color='r')
plt.axvline(280, linestyle = '--',color='g')
plt.axvline(0, linestyle = '--',color='g',label='1 bei 0 und 280')
plt.grid()
plt.ylim(-0.05,1.05)
plt.xlim(-10,290)
plt.legend(loc='best')
plt.savefig('Vorbereitungsplot1.pdf')

plt.clf()

plt.plot(L2, (1-(L2/roc)), label=r'$r_1 = 1400mm$ & $r_2=$eben')
plt.xlabel(r'Resonatorabstand $L/cm$')
plt.ylabel(r'Stabilitätsparameter $g_1 \cdot g_2$')
plt.axhline(1, linestyle = '--',color='r', label='Laserbereich')
plt.axhline(0, linestyle = '--',color='r')
plt.axvline(140, linestyle = '--',color='g')
plt.axvline(0, linestyle = '--',color='g',label='1 bei 0 und 140')
plt.grid()
plt.ylim(-0.05,1.05)
plt.xlim(-10,150)
plt.legend(loc='best')
plt.savefig('Vorbereitungsplot2.pdf')


