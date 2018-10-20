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
hr1 = np.mean(hr1)/60
hr2 = np.mean(hr2)/60
print("Calc Heating Rate1:",hr1)
print("Calced Heating Rate2:",hr2)
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
print('Parameter Fit 1.')
for p in params1:
    print(p)
#2
params , cov = curve_fit(fit , ut2 ,uc2 )
params2 = correlated_values(params, cov)
print('Parameter Fit 2.')
for p in params2:
    print(p)

## Prepare Plot Data ##
fc1 = np.append(mc1,current1)
fc1 = unp.uarray(fc1,np.absolute(fc1)*0.01)
fc2 = np.append(mc2,current2)
fc2 = unp.uarray(fc2,np.absolute(fc2*0.01))
xf1 = np.linspace(-90,55, 100)
print('vLines für Cut 1.:', ut1[uc1 == min(uc1)] ,ut1[ uc1 == min(uc1[ut1>-20])])
print('vLines für Cut 2.:', ut2[uc2 == min(uc2)] ,ut2[ uc2 == min(uc2[ut2>-20])])
#### Underground Plots ##
#plt.errorbar(temp1, noms(fc1), yerr = stds(fc1), fmt='xb', label = "Messdaten Heizrate = 2")
#plt.plot(ut1, uc1, 'or', label='Untergrund-Daten')
#plt.plot(xf1, fit(xf1,*noms(params1)), '--g', label='Untergrundausgleich')
#plt.vlines([-36.6,2.4],-11,31,colors='r',linestyle='--',label='Auswertungsbreich')
#plt.vlines(-10.2,-11,31,colors='k',linestyle='--',label='Tmax')
#plt.ylim(-10,30)
#plt.xlim(-75,55)
#plt.grid()
#plt.xlabel(r'$Probentemperatur \; / \;\degree C$')
#plt.ylabel(r'$Relaxationsstrom \; / \; A\cdot 10^{-11}$')
#plt.legend(loc='best')
##plt.show()
#plt.savefig('M1UGplot.pdf')
#plt.clf()
#plt.errorbar(temp2, noms(fc2), yerr = stds(fc2), fmt='xb', label = "Messdaten Heizrate = 1.5")
#plt.plot(ut2, uc2, 'or', label='Untergrund-Daten')
#plt.plot(xf1, fit(xf1,*noms(params2)), '--g', label='Untergrundausgleich')
#plt.vlines([-41.8,-0.4],-11,31,colors='r',linestyle='--',label='Auswertungsbreich')
#plt.vlines(-15.8,-11,31,colors='k',linestyle='--',label='Tmax')
#plt.ylim(-10,30)
#plt.xlim(-89,45)
#plt.grid()
#plt.xlabel(r'$Probentemperatur \; / \;\degree C$')
#plt.ylabel(r'$Relaxationsstrom \; / \; A\cdot 10^{-11}$')
#plt.legend(loc='best')
##plt.show()
#plt.savefig('M2UGplot.pdf')

## Chop Data thats not in the Space to Analyse & Subtraced the Underground ##
current1 = current1[(ctemp1 > ut1[uc1 == min(uc1)]) & (ctemp1 < ut1[uc1 == min(uc1[ut1>-20])])]
current2 = current2[(ctemp2 > ut2[uc2 == min(uc2)]) & (ctemp2 < ut2[ uc2 == min(uc2[ut2>-20])])]
ctemp1 = ctemp1[(ctemp1 > ut1[uc1 == min(uc1)]) & (ctemp1 < ut1[uc1 == min(uc1[ut1>-20])])]
ctemp2 = ctemp2[(ctemp2 > ut2[uc2 == min(uc2)]) & (ctemp2 < ut2[ uc2 == min(uc2[ut2>-20])])]
print('xpos Tmax 1:', ctemp1[current1==max(current1)])
print('xpos Tmax 2:', ctemp2[current2==max(current2)])
current1-= fit(ctemp1,*noms(params1))
current2-= fit(ctemp2,*noms(params2))

### Preps 4 Second Method for W Fit ### 
def nearest_toZero(array):
	idx = (np.abs(array - 0)).argmin()
	return array[idx]
print('1.i(T) ungefähr 0 bei  :', ctemp1[current1 == nearest_toZero(current1[ctemp1 > -10.2])])
print('2.i(T) ungefähr 0 bei  :', ctemp2[current2 == nearest_toZero(current2[ctemp2 > -15.8])])

Tstar1 = ctemp1[current1 == nearest_toZero(current1[ctemp1 > -10.2])]
Tstar2 = ctemp2[current2 == nearest_toZero(current2[ctemp2 > -15.8])]
def lnint(T,yt,hr):
	array = np.array([])
	for t in T:
		array = np.append(array,np.trapz(yt[T<=t],T[T<=t])/(yt[T == t]*hr))
	return array
T1 = ctemp1[ctemp1 < Tstar1] + const.zero_Celsius	
T2 = ctemp2[ctemp2 < Tstar2] + const.zero_Celsius	
yT1 = current1[ctemp1 < Tstar1]
yT2 = current2[ctemp2 < Tstar2]
#inti1 = (current1[ctemp1 < Tstar1],ctemp1[ctemp1 < Tstar1])
#inti2 = (current2[ctemp2 < Tstar2],ctemp2[ctemp2 < Tstar2])
#m2y1 = inti1/(current1*hr1)
#m2y2 = inti2/(current2*hr2)
m2y1 = lnint(T1,yT1,hr1)
m2y2 = lnint(T2,yT2,hr2)
m2x1 = T1[m2y1>0] 
m2x2 = T2[m2y2>0]
m2y1 = m2y1[m2y1>0]
m2y2 = m2y2[m2y2>0]
m2y1 = np.log(m2y1)
m2y2 = np.log(m2y2)
m2x1 = 1/m2x1
m2x2 = 1/m2x2
## First Method for W Fit ##
zeroCinK = const.zero_Celsius
invt1 = zeroCinK + ctemp1[ctemp1<ctemp1[current1==max(current1)]]
invt2 = zeroCinK + ctemp2[ctemp2<ctemp2[current2==max(current2)]] 
lnc1 = current1[ctemp1<ctemp1[current1==max(current1)]]
lnc2 = current2[ctemp2<ctemp2[current2==max(current2)]] 
invt1 = invt1[lnc1>0]
invt2 = invt2[lnc2>0]
m2c1i = lnc1
m2c2i = lnc2
lnc1 = np.log(lnc1[lnc1>0])
lnc2 = np.log(lnc2[lnc2>0])
invt1 = 1/invt1
invt2 = 1/invt2
def lnWfit(invt,A,B):
	return A*invt + B 
 
params , cov = curve_fit(lnWfit , invt1 ,lnc1)
paramslnW1 = correlated_values(params, cov)
print('Parameter lnWFit 1.')
for p in paramslnW1:
    print(p)

params , cov = curve_fit(lnWfit , invt2 ,lnc2)
paramslnW2 = correlated_values(params, cov)
print('Parameter lnWFit 2.')
for p in paramslnW2:
    print(p)
kbolz = const.k
print('1.Methode W mit Heizrate 2:', (-paramslnW1[0]/kbolz))
print('1.Methode W mit Heizrate 1,5 :', (-paramslnW2[0]/kbolz))

#### Plots of the First Method ###
#xt = np.linspace(0.0038,0.0043,100)
#plt.plot(invt1,lnc1, 'or', label='Ausgewertete Daten Heizrate = 2')
#plt.plot(invt2,lnc2, 'ob', label='Ausgewertete Daten Heizrate = 1.5')
#plt.plot(xt,lnWfit(xt,*noms(paramslnW1)),'--g', label='Ausgleichsgerade Heizrate = 2')
#plt.plot(xt,lnWfit(xt,*noms(paramslnW2)),'--k', label='Ausgleichsgerade Heizrate = 1.5')
#plt.grid()
##plt.ylim(-10,30)
##plt.xlim(,45)
#plt.xlabel(r'$Probentemperatur \; \frac{1}{T}\; / \;K$')
#plt.ylabel(r'$Relaxationsstrom  \; ln( I )\; / \; A\cdot 10^{-11}$')
#plt.legend(loc='best')
##plt.show()
#plt.savefig('1.MethFitW.pdf')

### Second Method ### 

params , cov = curve_fit(lnWfit , m2x1 ,m2y1)
params2lnW1 = correlated_values(params, cov)
print('Parameter Meth 2 lnWFit 1.')
for p in params2lnW1:
    print(p)

params , cov = curve_fit(lnWfit , m2x2 ,m2y2)
params2lnW2 = correlated_values(params, cov)
print('Parameter Meth 2 lnWFit 2.')
for p in params2lnW2:
    print(p)
print('2 Meth W Hr 2:', (-params2lnW1[0]/kbolz))
print('2 Meth W Hr 1.5 :',(-params2lnW2[0]/kbolz))
print('2 Meth teta Hr2:', (1/np.exp(noms(params2lnW1[1]))))
print('2 Meth teta Hr1.5 :', (1/np.exp(noms(params2lnW2[1]))))

### Plots 2 Methode ###
xt = np.linspace(0.0036,0.0044,100)
plt.plot(m2x1,m2y1, 'or', label='Daten Heizrate = 2')
plt.plot(m2x2,m2y2, 'ob', label='Daten Heizrate = 1.5')
plt.plot(xt,lnWfit(xt,*noms(params2lnW1)),'--g', label='Ausgleichsgerade Heizrate = 2')
#plt.plot(xt,lnWfit(xt,*noms(paramslnW2)),'--k', label='Ausgleichsgerade Heizrate = 1.5')
#plt.grid()
##plt.ylim(-10,30)
##plt.xlim(,45)
plt.xlabel(r'$Probentemperatur \; \frac{1}{T}\; / \;K$')
plt.ylabel(r'$Relaxationsstrom  \; ln ( \int_T ^{T*} I(T)dT\; / \; I(T)\cdot b ) \; $')
#plt.legend(loc='best')
plt.show()
#plt.savefig('1.MethFitW.pdf')

#Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',newline= r'\\'+'\n' )

