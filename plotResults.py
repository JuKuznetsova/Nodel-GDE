# plot calculation results
import math

import matplotlib.pyplot as plt
from math import log10


def plot_distribution(x, y):
    V_tot = 0.0
    N_tot = 0.0
    for i in range(len(x)):
        V_tot = V_tot + x[i] * y[i]
        N_tot = N_tot + y[i]
    print('V_tot = ', V_tot)
    print('N_tot = ', N_tot)

    for i in range(len(x)):
        y[i] = y[i] / N_tot
        x[i] = log10(x[i] * N_tot / V_tot)

    plt.plot(x, y)
    plt.savefig('fig_res_Coagulation.png')


# plot graphic x = nodel, y = PNC (Particle Number Concentration(#/m**3))
def plot_distribution_NodelToPNC(y):
    plt.plot(range(len(y)), y)
    plt.xlabel('Nodel')
    plt.ylabel('Particle Number Concentration(#/m**3)')
    plt.savefig('fig_res_CoagulationNucleation.png')


data = open('res_Coagulation_Nucleation.txt', 'r')
y = []
for line in data.readlines():
    y.append(float(line.split()[1]))
data.close()

plot_distribution_NodelToPNC(y)
print('Done!')
