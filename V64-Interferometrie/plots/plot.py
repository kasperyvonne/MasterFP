import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 100, 10**6)
kw, kmax , kmin = np.genfromtxt('Kontrast.txt',unpack = True) # grad , mV , V
k = (kmax[::] - kmin[::])/(kmax[::] + kmin[::])
umrechni = 2*np.pi/360
#def mittel(x):              #the real mean()-ing of life
#    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
#def relf(l,m):  #in Prozent
#    return (np.absolute(l-m)/l)*100
def fitf(x,a,b,c):
	return np.absolute(a*np.sin(b*x*umrechni)+ c)
#
#Fit
params , cov = curve_fit(fitf , kw ,k)
params = correlated_values(params, cov)
for p in params:
    print(p)
print(x[fitf(x,*noms(params)) == max(fitf(x,*noms(params)))])
plt.plot(kw,k, 'xr')
plt.plot(x,fitf(x,*noms(params)), 'g--')
plt.show()
#
#
##Tabelle
## np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )
##plt.subplot(1, 2, 1)
#plt.plot(x, y, label='Kurve')
#plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
#plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
#plt.legend(loc='best')
#plt.savefig('build/plot.pdf')
#plt.clf()
##plt.subplot(1, 2, 2)
#plt.plot(x, y, label='Kurve')
#plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
#plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
#plt.legend(loc='best')
#
## in matplotlibrc leider (noch) nicht möglich
##plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/plot2.pdf')
