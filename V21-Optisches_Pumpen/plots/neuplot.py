import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import correlated_values, correlation_matrix
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy import constants


def linear(x, m, b):
    return (x*m + b)


def Magnetfeld(I, N, R):
    return constants.mu_0 * (8*I*N)/(np.sqrt(125)*R)


bohr = const.value('Bohr magneton')
freq, sweep1, sweep2, hori1, hori2 = np.genfromtxt('Aufgabenteil_c.txt',
                                                   unpack='True')
verti = 0.1*2.31  # umedrehung *0.1V
freq = freq*1000  # in Hz
Isweep1 = 0.1*sweep1
Isweep2 = 0.1*sweep2

Ihori1 = 0.03*hori1
Ihori2 = 0.03*hori2

N_Sweep = 11
R_Sweep = 16.39*10**(-2)

N_Horizontal = 154
R_Horizontal = 15.79*10**(-2)

N_Vertical = 20
R_Vertical = 11.735*10**(-2)

B1_Sweep = Magnetfeld(Isweep1, N_Sweep, R_Sweep)
B1_Horizontal = Magnetfeld(Ihori1, N_Horizontal, R_Horizontal)
B1 = (B1_Sweep + B1_Horizontal)

B2_Sweep = Magnetfeld(Isweep2, N_Sweep, R_Sweep)
B2_Horizontal = Magnetfeld(Ihori2, N_Horizontal, R_Horizontal)
B2 = (B2_Sweep + B2_Horizontal)

paramsB1, covB1 = curve_fit(linear, freq, B1)
errorsB1 = np.sqrt(np.diag(covB1))
mB1 = ufloat(paramsB1[0], errorsB1[0])
bB1 = ufloat(paramsB1[1], errorsB1[1])

paramsB2, covB2 = curve_fit(linear, freq, B2)
errorsB2 = np.sqrt(np.diag(covB2))
mB2 = ufloat(paramsB2[0], errorsB2[0])
bB2 = ufloat(paramsB2[1], errorsB2[1])

g1 = constants.h/(mB1 * bohr)
g2 = constants.h/(mB2 * bohr)

print('Landésche gF')
print("g1: ", g1)
print("g2: ", g2)
g = g1/g2
print("Verhältnis 1/2: ", g)
print('')

print('---------------------------------------------------------------------')
print('Horizontalkomponenten')
print('Steigung 1: ', mB1)
print('Steigung 2: ', mB2)
print('Horizontalkomponente 1: ', bB1)
print('Horizontalkomponente 2: ', bB2)
print('')

J = 0.5
S = 0.5
L = 0
gJ = (3.0023 * (J**2 + J) + 1.0023 * ((S**2 + S)
      - (L**2 + L))) / (2 * (J**2 + J))
# I1 = gJ / (4 * g1) - 1 + unp.sqrt((gJ / (4 * g1) - 1)**2
#                                   + 3 * gJ / (4 * g1) - 3 / 4)
# I2 = gJ / (4 * g2) - 1 + unp.sqrt((gJ / (4 * g2) - 1)**2
# + 3 * gJ / (4 * g2) - 3 / 4)

I1 = 0.5*((gJ/g1)-1)
I2 = 0.5*((gJ/g2)-1)

print('--------------------------------------------------------------')
print('Kernspins')
print('Kernspin 1: ', I1)
print('Kernspin 2:', I2)
print('')

# assessment quadratic zeeman
U1 = g1*bohr*np.max(B1)+g1**2*bohr**2*np.max(B1)**2*(1-2*2)/(4.53e-24)
U2 = g2*bohr*np.max(B2)+g2**2*bohr**2*np.max(B2)**2*(1-2*3)/(2.01e-24)

print('------------------------------------------------------------')
print('Quadratische Zeeman-Aufspaltung')
print('Maximales BFeld1: ', np.round(np.max(B1)*10**6, 2))
print('Maximales BFeld2: ', np.round(np.max(B2)*10**6, 2))
print('Quadratische Zeeman-Aufspaltung 1 in eV: ', U1/constants.e)
print('Quadratische Zeeman-Aufspaltung 2 in eV: ', U2/constants.e)

x_plot = np.linspace(-100, 1050000)
plt.plot(freq, B1*10**6, 'bx', label='Isotop 1')
plt.plot(x_plot, linear(x_plot, *paramsB1)*10**6, 'b-', label='Ausgleichsgerade 1', linewidth=1)
plt.plot(freq, B2*10**6, 'rx', label='Isotop 2')
plt.plot(x_plot, linear(x_plot, *paramsB2)*10**6, 'r-', label='Ausgleichsgerade 2', linewidth=1)
plt.xlim(0, 1050000)
plt.xlabel(r'$f \:/\: $kHz')
plt.ylabel(r'$B \:/\: \mu}$T')
plt.legend(loc='best')
plt.savefig('neuBFelder.pdf')
plt.clf()
