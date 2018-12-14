import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
lambdvac = 632.990 #nano meter
lengL = ufloat(100,0.1) #milli meter
T = 22.8 #celsius
T+= const.zero_Celsius #kelvin
### Kontrast Kram
kw, kmax , kmin = np.genfromtxt('Kontrast.txt',unpack = True) # grad , mV , V
k = (kmax[::] - kmin[::])/(kmax[::] + kmin[::])
umrechni = 2*np.pi/360
def fitf(x,a,b,c):
	return np.absolute(a*np.sin(b*x*umrechni)+ c)
#
#Fit
params , cov = curve_fit(fitf , kw ,k)
params = correlated_values(params, cov)
for p in params:
    print(p)
x = np.linspace(0, 100, 10**6)
print(x[fitf(x,*noms(params)) == max(fitf(x,*noms(params)))])
plt.plot(kw,k, 'xr',label= "Kontrast")
plt.plot(x,fitf(x,*noms(params)), 'g--', label="Ausgleichskurve")
plt.xlabel(r'Polarisationswinkel $\theta$')
plt.ylabel(r'Kontrast $K$')
plt.legend(loc='best')
plt.grid()
plt.savefig('Kontrastfit.pdf')
plt.clf()
###Plätchen 
wd, c1,c2,c3,c4,c5,c6 = np.genfromtxt('plaettchen.txt', unpack = True)
def plaetbrech(tet, a):
	return a*tet #**2 
aes = []
ces = [c1,c2,c3,c4,c5,c6]

for c in ces: 
	params , cov = curve_fit(plaetbrech , wd ,c)
	params = correlated_values(params, cov)
	aes.append(params)
labels = ["Ausgleichskurve %i" % i for i in range(1,7) ]
for a,l in zip(aes,labels) :
	plt.plot(wd,plaetbrech(wd,noms(a)), label=l)
plt.plot(wd, c1,'o',alpha = 0.6 ,label='Messreihe 1')
plt.plot(wd, c2,'x',alpha = 0.6 ,label='Messreihe 2')
plt.plot(wd, c3,'v',alpha = 0.6 ,label='Messreihe 3')
plt.plot(wd, c4,'^',alpha = 0.6 ,label='Messreihe 4')
plt.plot(wd, c5,'<',alpha = 0.6 ,label='Messreihe 5')
plt.plot(wd, c6,'>',alpha = 0.6 ,label='Messreihe 6')
plt.xlabel(r'Winkel $\theta$')
plt.ylabel(r'Anzahl der Maxima $M$')
plt.legend(loc='best')
plt.grid()
plt.savefig('Plaettchenplot.pdf')
aes= unp.uarray(noms(aes),stds(aes))
brechisp = 1/(1-(lambdvac*aes /T))
np.savetxt('plaettab.txt',np.column_stack([noms(aes),stds(aes),noms(brechisp),stds(brechisp)]), delimiter=' & ',newline= r'\\'+'\n' )
### Luft 

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
