import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.signal import find_peaks
x = np.linspace(0, 10, 1000)
mhub = const.value('Bohr magneton') #das gelibete Borhsche Magneton zeigt wie man Scipy Constants benutzt
def mittel(x):              #the real mean()-ing of life
    return ufloat(np.mean(x),np.std(x,ddof=1)/np.sqrt(len(x)))
def relf(l,m):  #in Prozent
    return (np.absolute(l-m)/l)*100
##Fit
#params , cov = curve_fit(f , x ,y )
#params = correlated_values(params, cov)
#for p in params:
#    print(p)
### T1 Bestimmung 
tau1, U1 = np.genfromtxt('T1Messung.txt', unpack = True)
U1*= 10**(-3) # umrechnen in Volt
def T1fit(tau,m0, t1,b):
	return m0*(1-2*np.exp(-(tau/t1)))+b 
params , cov = curve_fit(T1fit , tau1 ,U1)
params = correlated_values(params, cov)
print("T1 fit:")
pams=['U0:','T1:','b:']
for i,j in zip(params,pams):
	print(j,i)
T1 = params[1]
taufit1 = np.linspace(min(tau1-1),max(tau1+1),1000)
plt.plot(tau1, U1,'rx', label='Induzierte Spannung $U_z$')
plt.plot(taufit1, T1fit(taufit1,*noms(params)), 'b--',label='Ausgleichskurve')
plt.axhline(y=0, xmin =min(tau1-1),xmax =max(tau1+1), color = 'k')
plt.axhline(y=noms(params[0]), xmin =min(tau1-1),xmax =max(tau1+1), color = 'g', linestyle ='--', label = '$U_0$')
plt.axhline(y=-noms(params[0]), xmin =min(tau1-1),xmax =max(tau1+1), color = 'g', linestyle ='--')
#plt.xlim(min(tau1-0.5),max(tau1+0.5))
plt.xlabel(r'$\tau /$s')
plt.ylabel(r'$U_z /$V')
plt.legend(loc='best')
plt.xscale("log")
plt.grid()
plt.savefig('T1plot.pdf')
#plt.show()
plt.clf()

### Meiboom-Gill-Methode
tmg, Umg = np.genfromtxt('MGM.txt', unpack = True)
xmgpeaks, ymgpeaks = find_peaks(Umg)
fUmg = Umg[xmgpeaks]
ftmg = tmg[xmgpeaks]
ftmg = ftmg[fUmg > 306]
fUmg = fUmg[fUmg > 306]

def T2fit(t, lnu0, t2 ):
	return lnu0-(t/t2)
params , cov = curve_fit(T2fit , ftmg ,np.log(fUmg))
params = correlated_values(params, cov)
print("T2 fit:")
pams=['lnU0:','T2:' ]
for i,j in zip(params,pams):
	print(j,i)
print('U0:', np.exp(noms(params[0])))
T2 = params[1]

#plt.plot(tmg[xmgpeaks], Umg[xmgpeaks], 'rx')
tmgfit = np.linspace(-100,1100, 1000)
plt.plot(ftmg,np.log(fUmg), 'rx', label = 'logarith. Signalspitzen')
plt.plot(tmgfit,T2fit(tmgfit,*noms(params)), 'b--', label = 'Ausgleichsgerade')
plt.xlim( xmin =min(ftmg-50),xmax =max(ftmg+50))
plt.xlabel(r'$t/$ms')
plt.ylabel(r'log$(U_y /$V)')
plt.legend(loc='best')
plt.grid()
plt.savefig('T2plot.pdf')
#plt.show()
plt.clf()

### Diffusionskonstante
dtau, dU, thalb = np.genfromtxt('DiffMessung.txt', unpack = True)
dtau*=2*10**(-3) 
dU*= -10**(-3)
thalb*=10**(-6)
mth = mittel(thalb)
print('### thalbe:',mth)
d = 4.4 # mm Proben Durchmesser
d*=10**(-3)
A = (4*2.2)/(d*mth)
gamma = ufloat(425.77469 , 0.000013)*10**6 # gamma/2pi = Hz/Tesla
print('Gamma in Einheiten von 2pi:', gamma)
print('Gradient G in Einheiten von 2pi', A/gamma)  
def Dfit(t,u0,D):
	return u0-t/noms(T2)-(1/12)*D*(noms(A)**2)*(t**3)	
#def Dfit(t,u0,D, u1):
#	return  u0*np.exp(-(1/12)*D*(noms(A)**2)*(t**3)-t/noms(T2)) +u1
params , cov = curve_fit(Dfit ,dtau,np.log(dU ))
params = correlated_values(params, cov)
print("D fit:")
pams=['lnU0:', 'D:']
for i,j in zip(params,pams):
	print(j,i)
Diff = params[1]
print('U0:', np.exp(noms(params[0])))
dtaufit = np.linspace(0,0.035,1000)
plt.plot(dtau, dU, 'rx', label ='logarith. Spin Echo')
plt.plot(dtaufit,Dfit(dtaufit,*noms(params)),'b--', label = 'Ausgleichskurve')
plt.xlim(0,0.035)
plt.xlabel(r'$t = 2 \tau/$s')
plt.ylabel(r'log$(U_y /$V)')
plt.legend(loc='best')
plt.grid()
plt.savefig('Dplot.pdf')
#plt.show()
plt.clf()
###Viskosi Messung
vt,vd = np.genfromtxt('Visdata.txt', unpack=True)
def visfit(x,a,b):
	return a*x +b
vist = 15.29*60
print('vis t =', vist)
params , cov = curve_fit(visfit ,vt,vd)
params = correlated_values(params, cov)
delta = visfit(vist,*params)
rho = 1000
alpha = 1.024*10**(-9)
vis = rho * alpha *(vist - delta)
print('delta vis :',delta)
print('viskosität:', vis)
###Berechnung des Molekuelradiuses:
rstokes = (const.k * (const.zero_Celsius + 20))/(6* np.pi * vis * Diff)
rhcp = ((0.74*3*const.value('atomic mass constant')*(16+2*1))/(4*rho*np.pi))**(1/3) 
rkrit = ((3*const.k*647.05)/(128*np.pi*22.04*10**6))**(1/3)
print(const.R)
print('kritischer Druck 22.04 MPa und krit Temp 647.05 Kelvin')
print('Stokes Molrad:', rstokes)
print('hcp molrad:', rhcp)
print('krit molrad:',rkrit)
######Diskussion:
print('######Diskussion:')
print('lit für T1=T2=2500ms')
print('relf T1',relf(2500,T1*10**(3)))
print('relf T2',relf(2500,T2))
print('Vgl rstokes mit hcp und krit:')
print('stokes hcp:',relf(rhcp,rstokes))
print('stokes krit:', relf(rkrit,rstokes))
print('lit Wert D:',2.023*10**(-9))
print('relf D:', relf(2.023*10**(-9),Diff))
print('lit Wert vis:', 0.01*10**(-1))
print('relf vis:', relf(0.01*10**(-1),vis))
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
