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
    plt.savefig('fig_res_CoagulationDistribution.png')


# x = []
# y = []
# data = open('res_Coagulation_0.01s.txt', 'r')
# i = 0
# for line in data.readlines():
#     if i % 2 == 0:
#         x.append(float(line.split()[0]))
#         y.append(float(line.split()[1]))
#     i = i+1
# data.close()
#
# plot_distribution(x, y)
# print('Done!')
