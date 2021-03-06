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


def Helmhotz(N, I, R):
    return ((const.mu_0*(8*N*I))/(np.sqrt(125)*R))


def g(x, m, b):
    return (m*x+b)


def h(x, a, b, c):
    return a + b/(x+c)


bohr = const.value('Bohr magneton')
freq, sweep1, sweep2, hori1, hori2 = np.genfromtxt('Aufgabenteil_c.txt',
                                               unpack='True')
Vpp1, Peaks1, delT1 = np.genfromtxt('Aufgabenteil_i1.txt', unpack='True')
Vpp2, Peaks2, delT2 = np.genfromtxt('Aufgabentei_i2.txt', unpack='True')

verti = 0.1*2.31  # umedrehung *0.1V
sweep1 = 0.1*sweep1
sweep2 = 0.1*sweep2

hori1 = 0.03*hori1
hori2 = 0.03*hori2

#print(sweep1, sweep2, hori1, hori2)
Nhori = 154
Rhori = 0.1579  # m

Nsweep = 11
Rsweep = 0.1639  # m

Nverti = 20
Rverti = 0.11735

Bsweep1 = Helmhotz(Nsweep, sweep1, Rsweep)
Bsweep2 = Helmhotz(Nsweep, sweep2, Rsweep)
Bhori1 = Helmhotz(Nhori, hori1, Rhori)
Bhori2 = Helmhotz(Nhori, hori2, Rhori)
# print(Bsweep1*10**6)
# print(Bsweep2*10**6)
# print(Bhori1*10**6)
# print(Bhori2*10**6)

B1 = (Bsweep1 + Bhori1)
B2 = (Bsweep2 + Bhori2)

# print(B1)
# print(B2)
B2fit = B2[0:]
freqfit = freq[0:]
x_plot = np.linspace(-100, 1050)
params1, covariance1 = curve_fit(g, freq, B1)
errors1 = np.sqrt(np.diag(covariance1))

params2, covariance2 = curve_fit(g, freq, B2)
errors2 = np.sqrt(np.diag(covariance2))

m1 = ufloat(params1[0], errors1[0])
m2 = ufloat(params2[0], errors2[0])

g1 = const.h/(m1*bohr)*1000
g2 = const.h/(m2*bohr)*1000

Verhaeltnisg = g1/g2

J = 0.5
S = 0.5
L = 0
gJ = (3.0023 * (J**2 + J) + 1.0023 * ((S**2 + S)
      - (L**2 + L))) / (2 * (J**2 + J))
#I1 = gJ / (4 * g1) - 1 + unp.sqrt((gJ / (4 * g1) - 1)**2
#                                  + 3 * gJ / (4 * g1) - 3 / 4)
#I2 = gJ / (4 * g2) - 1 + unp.sqrt((gJ / (4 * g2) - 1)**2
#                                  + 3 * gJ / (4 * g2) - 3 / 4)
I1 = 0.5*((gJ/g1)-1)
I2 = 0.5*((gJ/g2)-1)

plt.plot(freq, B1*10**6, 'bx', label='Isotop 1')
plt.plot(x_plot, g(x_plot, *params1)*10**6, 'b-', label='Ausgleichsgerade 1', linewidth=1)
plt.plot(freq, B2*10**6, 'rx', label='Isotop 2')
plt.plot(x_plot, g(x_plot, *params2)*10**6, 'r-', label='Ausgleichsgerade 2', linewidth=1)
plt.xlim(0, 1050)
plt.xlabel(r'$f \:/\: $kHz')
plt.ylabel(r'$B \:/\: \mu}$T')
plt.legend(loc='best')
plt.savefig('BFelder.pdf')
plt.clf()

# quadratischer Zeemaneffekt

U1 = g1*bohr*B1[-1]+g1**2*bohr**2*B1[-1]**2*(1-2*2)/(4.53e-24)
U2 = g2*bohr*B2[-1]+g2**2*bohr**2*B2[-1]**2*(1-2*3)/(2.01e-24)

# letzter Teil
T1 = Peaks1 / delT1
T2 = Peaks2 / delT2

params3, covariance3 = curve_fit(h, Vpp1, T1)
errors3 = np.sqrt(np.diag(covariance3))

params4, covariance4 = curve_fit(h, Vpp2, T2)
errors4 = np.sqrt(np.diag(covariance4))

bt1 = ufloat(params3[1], errors3[1])
bt2 = ufloat(params4[1], errors4[1])
bverh = bt1/bt2
print(bt1)
print(bt2)
xtrans = np.linspace(0.5, 10.5)
plt.plot(Vpp1, T1, 'x', label=r'$\mathrm{^{87}Rb}$')
plt.plot(xtrans, h(xtrans, *params3), 'b-', label='Ausgleichskurve', linewidth=1)
plt.xlabel(r'$Amplitude \:/\: $V')
plt.ylabel(r'$T \:/\:}$ms')
plt.legend(loc='best')
plt.savefig('Trans1.pdf')
plt.clf()
plt.plot(Vpp2, T2, 'x', label=r'$\mathrm{^{85}Rb}$')
plt.plot(xtrans, h(xtrans, *params4), 'b-', label='Ausgleichskurve', linewidth=1)
plt.xlabel(r'$Amplitude \:/\: $V')
plt.ylabel(r'$T \:/\:}$ms')
plt.legend(loc='best')
plt.savefig('Trans2.pdf')
plt.clf()

Ampli1 = ufloat(70.51, 0.7)
Ampli2 = ufloat(149.33, 1.5)
print("Das angelegte vertikale Feld entspricht:",
      Helmhotz(Rverti, Nverti, verti), "Tesla.")
print("Die Fitparameter für Isotop1 sind m=", params1[0], '±', errors1[0],
      "und b=", params1[1], '±', errors1[1])
print("Die Fitparameter für Isotop2 sind m=", params2[0], '±', errors2[0],
      "und b=", params2[1], '±', errors2[1],)
print("Landefaktor1=", g1)
print("Landefaktor2=", g2)
print("Das Verhältnis der beiden Faktoren ist", Verhaeltnisg)
print(gJ)
print("Der Kernsprin für Isotop1 ist:", I1)
print("Der Kernsprin für Isotop2 ist:", I2)
print("quadratische I1:", U1)
print("quadratische I2:", U2)
print("b Verhältnis (soll: 1,5)", bverh)
print(Ampli1/Ampli2)
print(B1[-1])
print(B2[-1])
print("params transfit1:")
print(ufloat(params3[0], errors3[0]))
print(ufloat(params3[1], errors3[1]))
print(ufloat(params3[2], errors3[2]))
print("params transfit2:")
print(ufloat(params4[0], errors4[0]))
print(ufloat(params4[1], errors4[1]))
print(ufloat(params4[2], errors4[2]))
print((ufloat(params3[1], errors3[1]))/(ufloat(params4[1], errors4[1])))
print("relative Abweichung ver.Magnetfeld:")
print((45.11-35.4)/45.11)
print("relative Abweichung Kernspin:")
print((I1-1.5)/1.5)
print((I2-2.5)/2.5)
print((ufloat(0.47,0.01)-ufloat(0.368, 0))/ufloat(0.368, 0))
print((ufloat(1.2,0.3)-ufloat(1.5, 0))/ufloat(1.5, 0))



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
