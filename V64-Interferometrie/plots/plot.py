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
lambdvac = 632.990*10**(-9) #nano meter
lengL = ufloat(100,0.1)*10**(-3) #milli meter
T = 22.8 #celsius
T+= const.zero_Celsius #kelvin
d = 1*10**(-3)
### Kontrast Kram
kw, kmax , kmin = np.genfromtxt('Kontrast.txt',unpack = True) # grad , mV , V
k = (kmax[::] - kmin[::])/(kmax[::] + kmin[::])
umrechni = 2*np.pi/360
def fitf(x,a,b,c):
	return np.absolute(a*np.sin(2*x*umrechni + b)+ c)
#
#Fit
np.savetxt('kontrasttab.txt',np.column_stack([kw,noms(k)]), delimiter=' & ',newline= r'\\'+'\n' )
params , cov = curve_fit(fitf , kw ,k, [0.6 , 0.1 , 0])
params = correlated_values(params, cov)
for p in params:
    print(p)
x = np.linspace(0, 100, 10**6)
print(x[fitf(x,*noms(params)) == max(fitf(x,*noms(params)))])
plt.plot(kw,k, 'xr',label= "Kontrast")
plt.plot(x,fitf(x,*noms(params)), 'g--', label="Ausgleichskurve")
plt.xlabel(r'Polarisationswinkel $\theta°$')
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
wd*= 2*np.pi/360		#umrechnen in Bogenmaß
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
plt.xlabel(r'Winkeländerung $\delta\theta /$rad')
plt.ylabel(r'Anzahl der Maxima $M$')
plt.legend(loc='best')
plt.grid()
plt.savefig('Plaettchenplot.pdf')
plt.clf()
aes= unp.uarray(noms(aes),stds(aes))
brechisp = 1/(1-(lambdvac*aes /(2*(20*np.pi/360)*d)))
print('Plaettchen Mittelbrech:', np.mean(brechisp))
np.savetxt('plaettab.txt',np.column_stack([noms(aes),stds(aes),noms(brechisp),stds(brechisp)]), delimiter=' & ',newline= r'\\'+'\n' )
### Luft 
druck, l1,l2,l3,l4 = np.genfromtxt('Luft.txt', unpack = True)
druck*=100 # mbar zu Pa
laes = [l1,l2,l3,l4]
brechis = [ ((l*lambdvac)/lengL) +1 for l in laes] #*10^9
def lorlorlaw(p,a): #taylored
	return (1+ a*p)
bpams=[]
for b in brechis:
	params , cov = curve_fit(lorlorlaw , druck ,noms(b)**2)
	params = correlated_values(params, cov)
	bpams.append(params)
llabels = ["Ausgleichskurve %i" % i for i in range(1,5) ]
d = np.linspace(8000,110000,10**3)
for a,l in zip(bpams,llabels):
	plt.plot(d,lorlorlaw(d,noms(a)), label=l)
l2labels = ["Messreihe %i" % i for i in range(1,5) ]
for b,l in zip(brechis,l2labels):
	plt.errorbar(druck,noms(b)**2,yerr = stds(b),fmt='x', label=l)
plt.xlabel(r'Druck $p /$Pa')
plt.ylabel(r'Brechungsindex $n^2$')
plt.xlim(min(druck)-10**3,max(druck)+10**3)
plt.legend(loc='best')
plt.grid()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.show()
plt.savefig('Luftplot.pdf')
np.savetxt('lufttab.txt',np.column_stack([noms(bpams),stds(bpams)]), delimiter=' & ',newline= r'\\'+'\n' )
print(np.mean(bpams))
print(((np.mean(bpams)*T*1013*10**2/(15+const.zero_Celsius)) + 1)**0.5)
print('litW:',(27663.8*10**(-8) +1))
print('relf Lit Luft:', relf((27663.8*10**(-8) ),np.sqrt(noms((np.mean(bpams)*T*1013*10**2/(15+const.zero_Celsius)) + 1)) -1))
print('relf Lit Plaettchen:', relf(.5,np.mean(brechisp)-1))



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
