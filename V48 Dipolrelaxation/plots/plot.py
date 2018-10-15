import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
x = np.linspace(0, 10, 1000)
mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
### Read in Data ###
temp1, current1 = np.genfromtxt('M1_Heizrate2.txt', unpack = True) #Messreihe 1 Heizrate 2
temp2, current2 = np.genfromtxt('M2_Heizrate15.txt', unpack = True) #Messreihe 2 Heizrate 1.5
#temp in celsius
#current in 10^-11 Ampere
### Flip Current ###
current1*=-1
current2*=-1 
### calc Heating Rate ###
hr1 = list(map(lambda x,y: x-y, temp1[1:], temp1[:-1]))
hr2 = list(map(lambda x,y: x-y, temp2[1:], temp2[:-1]))
hr1 = np.mean(hr1)
hr2 = np.mean(hr2)
print("Calc Heating Rate:",hr1)
print("Calced Heating Rate:",hr2)
### substract Offset ###
offs = 1.2 #*10^-11 Ampere
current1-= offs
current2-= offs
### Chop off messy Data ###
mc1 = current1[temp1<-50] #messy Data 1
mc2 = current2[temp2<= -50] #messy Data 2
current1 = current1[temp1>-50]
current2 = current2[temp2>-50]
ctemp1 = temp1[temp1 > -50]
ctemp2 = temp2[temp2 > -50]
### Make Underground Data ###
uc1 = np.array([])
ut1 = np.array([])
filwid = 5
for index, item in enumerate(current1):
	if current1[index + filwid] == current1[-1]:
		break
	if index < filwid:
		continue
	if (item < current1[index +filwid]) and (item < current1[index-filwid]):
		uc1 = np.append(uc1,item)
		ut1 = np.append(ut1,ctemp1[index])
filwod= 5
uc2 = np.array([])
ut2 = np.array([])
for index, item in enumerate(current2):
	if current2[index + filwod] == current2[-1]:
		break
	if index < filwod:
		continue
	if (item < current2[index +filwod]) and (item < current2[index-filwod]):
		uc2 = np.append(uc2,item)
		ut2 = np.append(ut2,ctemp2[index])
### Underground Fit ###
def fit(T,A,B,C):
	return A*np.exp(C*T) + B
#1

params , cov = curve_fit(fit , ut1 ,uc1 )
params1 = correlated_values(params, cov)
for p in params1:
    print(p)
#2
params , cov = curve_fit(fit , ut2 ,uc2 )
params2 = correlated_values(params, cov)
for p in params2:
    print(p)

## Prepare Plot Data ##
fc1 = np.append(mc1,current1)
fc1 = unp.uarray(fc1,np.absolute(fc1)*0.01)
fc2 = np.append(mc2,current2)
fc2 = unp.uarray(fc2,np.absolute(fc2*0.01))
xf1 = np.linspace(-90,55, 100)

#plt.errorbar(temp1, noms(fc1), yerr = stds(fc1), fmt='xb', label = "Messdaten Heizrate = 2")
#plt.plot(ut1, uc1, 'or', label='Untergrund-Daten')
#plt.plot(xf1, fit(xf1,*noms(params1)), '--g', label='Untergrundausgleich')
#plt.ylim(-10,30)
#plt.xlabel(r'$Probentemperatur \; / \;\degree C$')
#plt.ylabel(r'$Relaxationsstrom \; / \; A\cdot 10^{-11}$')
#plt.legend(loc='best')
#plt.show()
#plt.savefig('.pdf')
#plt.clf()
plt.errorbar(temp2, noms(fc2), yerr = stds(fc2), fmt='xb', label = "Messdaten Heizrate = 1.5")
plt.plot(ut2, uc2, 'or', label='Untergrund-Daten')
plt.plot(xf1, fit(xf1,*noms(params2)), '--g', label='Untergrundausgleich')
plt.ylim(-10,30)
plt.xlabel(r'$Probentemperatur \; / \;\degree C$')
plt.ylabel(r'$Relaxationsstrom \; / \; A\cdot 10^{-11}$')
plt.legend(loc='best')
plt.show()
#plt.savefig('build/plot2.pdf')
#
#Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )

