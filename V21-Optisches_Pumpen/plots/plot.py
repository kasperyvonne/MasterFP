import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)


def mittel(x):              # the real mean()-ing of life
    return ufloat(np.mean(x), np.std(x, ddof=1)/np.sqrt(len(x)))


def relf(l, m):  # in Prozent
    return (np.absolute(l-m)/l)*100


def Helmhotz(R, N, I):
    return (const.mu_0*(8/np.sqrt(125))*(N/R)*I)


def Gerade(m, x, b):
    return (m*x+b)


x = np.linspace(0, 10, 1000)
mhub = const.value('Bohr magneton')
freq, sweep1, sweep2, hori1, hori2 = np.genfromtxt('Aufgabenteil_c.txt',
                                                   unpack='True')
verti = 0.1*2.31
sweep1 = 0.1*sweep1
sweep2 = 0.1*sweep2

hori1 = 0.3*hori1
hori2 = 0.3*hori2

Nhori = 154
Rhori = 0.1579  # m

Nsweep = 11
Rsweep = 0.1639  # m

Nverti = 20
Rverti = 0.11735

Bsweep1 = Helmhotz(Rsweep, Nsweep, sweep1)
Bsweep2 = Helmhotz(Rsweep, Nsweep, sweep2)
Bhori1 = Helmhotz(Rsweep, Nsweep, hori1)
Bhori2 = Helmhotz(Rsweep, Nsweep, hori2)

B1 = Bsweep1 + Bhori1
B2 = Bsweep2 + Bhori2
params1 , cov1 = curve_fit(gerade , x ,y )
params1 = correlated_values(params1, cov1)

plt.plot(freq, B1, 'x', label='Isotop 1')
plt.plot(freq, B2, 'x', label='Isotop 2')
plt.xlabel(r'$f \:/\: $kHz')
plt.ylabel(r'$B \:/\: \mu}$T')
plt.legend(loc='best')
plt.savefig('Testplot1.pdf')
plt.clf()

print("Das angelegte vertikale Feld entspricht:",
      Helmhotz(Rverti, Nverti, verti))

# Fit
# params , cov = curve_fit(f , x ,y )
# params = correlated_values(params, cov)
# for p in params:
#     print(p)
#
#
# #Tabelle
# np.savetxt('tab.txt',np.column_stack([x,y]), delimiter=' & ',
#             newline= r'\\'+'\n' )
# #plt.subplot(1, 2, 1)
# plt.plot(x, y, label='Kurve')
# plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
# plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
# plt.legend(loc='best')
# plt.savefig('build/plot.pdf')
# plt.clf()
# #plt.subplot(1, 2, 2)
# plt.plot(x, y, label='Kurve')
# plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
# plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
# plt.legend(loc='best')
#
# # in matplotlibrc leider (noch) nicht möglich
# #plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/plot2.pdf')
