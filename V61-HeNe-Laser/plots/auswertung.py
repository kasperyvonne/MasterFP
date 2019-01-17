import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
# import uncertainties.unumpy as unp
# import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy import constants
from scipy.signal import find_peaks

# Use latex fonts and text
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

# pre-assignment

#
# def bed1(L, r1, r2):
#     return ((1 - L / r1) * (1 - L / r2))
#
#
# def bed2(L, r1):
#     return (1 - L / r1)
#
#
# r1_1 = 1.4
# r2_1 = 1.4
#
# r1_2 = 1.4
#
# x = np.linspace(0, 3)
# a = bed1(x, r1_1, r2_1)
# b = bed2(x, r1_2)
#
# plt.figure(1)
# plt.plot(x, bed1(x, r1_1, r2_1), 'r-', label=r"$r_1$ = 1,4 m, $r_2$ = 1,4 m")
# plt.plot(x, bed2(x, r1_2), 'b-', label=r"$r_1$ = 1,4 m, $r_2$ = $\infty$ m")
# plt.axhline(y=1, color='k', linestyle='--')
# plt.axhline(y=0, color='k', linestyle='--')
# plt.xlabel(r'Resonatorlaenge $L$ / m')
# plt.ylabel(r'Stabilitaetsparameter $g_1 g_2$')
# plt.ylim(-1, 1.7)
# plt.xlim(0, 2.8)
# plt.legend(loc="best")
# plt.tight_layout()
# plt.savefig('Auswertung/Plots/g1g2.pdf')
# plt.clf()
#
# print('----------------------------------------------------------------------')

# Wert in m
l1 = 0.541
l2 = 1.055
l3 = 1.905

# Werte in MHz
P1 = np.array([75, 199, 274, 353, 473, 548, 623, 746])
P2 = np.array([143, 281, 424, 566, 705, 848, 986, 1129])
P3 = np.array([79, 158, 236, 315, 394, 473, 551, 626, 705, 784])
#--------Unsere Werte bis hier------
P1 = P1 * 10**6
P2 = P2 * 10**6
P3 = P3 * 10**6

nu1 = np.empty(len(P1) - 1)
nu2 = np.empty(len(P2) - 1)
nu3 = np.empty(len(P3) - 1)

for i in range(len(P1) - 1):
    nu1[i] = P1[i + 1] - P1[i]
for i in range(len(P2) - 1):
    nu2[i] = P2[i + 1] - P2[i]
for i in range(len(P3) - 1):
    nu3[i] = P3[i + 1] - P3[i]

dnu1 = ufloat(np.mean(nu1), sem(nu1))
dnu2 = ufloat(np.mean(nu2), sem(nu2))
dnu3 = ufloat(np.mean(nu3), sem(nu3))

# Some constants and values in SI units
c = constants.c
k = constants.k
T = 300
M = 20.1797 * constants.physical_constants["atomic mass constant"][0]
lam = 632.8 * 10**(-9)
f = c / lam

# v_mean from boltzman formular
v_mean = np.sqrt(2 * k * T / M)

print(v_mean)

# delta factor from doppler-effect
delta = (c + v_mean) / (c - v_mean)

# print frequency from 632.8nm wavelength and resulting doppler-effect
print("f_632,8 nm, delta_f_Doppler, lambda_mit_Doppler")
print(f, np.abs(2 * (f - delta * f)), np.abs(c / (delta * f)) * 10**9)

# print mean frequency between modes
print("delta_nu1, delta_nu2, delta_nu3:")
print(dnu1, dnu2, dnu3)

print('----------------------------------------------------------------------')
IDunkel = 1


def M00_Intensitaet(x, I0, sigma, x0):
    return I0 * np.exp(-2 * (x - x0)**2 / sigma**2)


def M01_Intensitaet(x, I0, sigma, x0):
    return  I0*(x -x0)**2 * np.exp(-2 * (x - x0)**2 / sigma**2)


x00, M00 = np.genfromtxt('TEM00.txt', unpack=True)
x01, M01 = np.genfromtxt('TEM01.txt', unpack=True)

# np.savetxt("M00_Tabelle.tex",
#            np.column_stack([
#                            x00[:30],
#                            M00[:30],
#                            x00[30:60],
#                            M00[30:60]
#                            ]), delimiter=' & ', newline=r' \\' + '\n',
#            fmt='%.1f & %.1f & %.1f & %.1f')
#
# np.savetxt("M01_Tabelle.tex",
#            np.column_stack([
#                            x01[:24],
#                            M01[:24],
#                            x01[24:],
#                            M01[24:]
#                            ]), delimiter=' & ', newline=r' \\' + '\n',
#            fmt='%.1f & %.1f & %.1f & %.1f')

M01 = M01
M00 = M00

params00, cov00 = curve_fit(M00_Intensitaet, x00, M00)
errors00 = np.sqrt(np.diag(cov00))
I00 = ufloat(params00[0], errors00[0])
sigma00 = ufloat(params00[1], errors00[1])
x0_00 = ufloat(params00[2], errors00[2])

params01, cov01 = curve_fit(M01_Intensitaet, x01, M01,  p0=(300, 1, 1))
errors01 = np.sqrt(np.diag(cov01))
I01 = ufloat(params01[0], errors01[0])
sigma01 = ufloat(params01[1], errors01[1])
x0_01 = ufloat(params01[2], errors01[2])

print("Mode: I0, sigma, x0")
print("TEM_00: ", I00, sigma00, x0_00)
print("TEM_01: ", I01, sigma01, x0_01)

xspace = np.linspace(min(x00) - 2.5, max(x00) + 2.5, 10000)

plt.figure(2)
plt.plot(x00, M00, 'rx', label=r"Messwerte TEM$_{00}$-Mode")
plt.plot(xspace, M00_Intensitaet(xspace, *params00), 'r-', label="Fit")
plt.xlabel(r'Position / mm')
plt.ylabel(r'Intensitaet / nA')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('M00.pdf')
plt.clf()

plt.figure(3)
plt.plot(x01, M01, 'rx', label=r"Messwerte TEM$_{01}$-Mode")
plt.plot(xspace, M01_Intensitaet(xspace, *params01), 'r-', label="Fit")
#plt.axhline(350, linestyle = '--',color='r')
plt.xlabel(r'Position / mm')
plt.ylabel(r'Intensitaet / nA')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('M01.pdf')
plt.clf()

print('----------------------------------------------------------------------')

phi, I = np.genfromtxt('Polarisation.txt', unpack=True)

# np.savetxt("Auswertung/Tabellen/Winkel_Tabelle.tex",
#            np.column_stack([
#                            phi,
#                            I
#                            ]), delimiter=' & ', newline=r' \\' + '\n',
#            fmt='%.0f & %.2f')


def Winkel_I(phi, A0, phi0):
    return A0 * np.cos((phi + phi0) * (np.pi / 180))**2

paramsWinkel, covWinkel = curve_fit(Winkel_I, phi, I)
errorsWinkel = np.sqrt(np.diag(covWinkel))
A0_Winkel = ufloat(paramsWinkel[0], errorsWinkel[0])
phi0_Winkel = ufloat(paramsWinkel[1], errorsWinkel[1])
# c_Winkel = ufloat(paramsWinkel[2], errorsWinkel[2])

phi_space = np.linspace(0, 360, 1000)

print("A0, phi0, c")
print(A0_Winkel, phi0_Winkel+360)

plt.figure(3)
plt.plot(phi, I, 'rx', label=r"Messwerte Polarisationsmessung")
plt.plot(phi_space, Winkel_I(phi_space, *paramsWinkel), 'r-', label="Fit")
plt.xlabel(r'Winkel / Grad')
plt.ylabel(r'Intensitaet / $\mu$A')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Polarisation.pdf')
plt.clf()

print('----------------------------------------------------------------------')
g = 1 / 100
#IDunkel = 0.00518
l = 53

x, I = np.genfromtxt('Wellenlängen.txt', unpack=True)

# np.savetxt("Auswertung/Tabellen/WL_Tabelle.tex",
#            np.column_stack([
#                            x[:22],
#                            I[:22],
#                            x[22:],
#                            I[22:]
#                            ]), delimiter=' & ', newline=r' \\' + '\n',
#            fmt='%.1f & %.3f & %.1f & %.3f')

#I = I - IDunkel

Peaks, _ = find_peaks(x=I)
Peak_Pos = x[Peaks]
Peaks = I[Peaks]

print("Peakpositionen: ", Peak_Pos)

plt.figure(4)
plt.plot(x, I, 'rx', label=r"Messwerte Beugungsbild")
plt.plot(Peak_Pos, Peaks, 'bx', label='Maxima')
# plt.plot(phi_space, Winkel_I(phi_space, *paramsWinkel), 'r-', label="Fit")
plt.xlabel(r'Winkel / Grad')
plt.ylabel(r'Intensitaet / $\mu$A')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Wellenlaenge.pdf')
plt.clf()

d1 = np.abs(Peak_Pos[3] - Peak_Pos[2])
d2 = np.abs(Peak_Pos[3] - Peak_Pos[4])

print("Abstände zum 0. Maximum: ", d1, "0", d2)

print("Winkel: ", np.arctan(d1 / l), 0, np.arctan(d2 / l))

lam1 = g * np.sin(np.arctan(d1 / l)) * 10**6
lam2 = g * np.sin(np.arctan(d2 / l)) * 10**6

print("Wellenlänge links, Wellenlänge rechts")
print(lam1, lam2)
print((lam1+lam2)/2)
print('----------------------------------------------------------------------')

print('Alles ausgeführt!')

print('----------------------------------------------------------------------')
