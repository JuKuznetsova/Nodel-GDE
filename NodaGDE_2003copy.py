"""
Эта программа соответствует коду, предложенному  и реализованному на языке С в статье
A Simple Numerical Algorithm and Software for Solution of Nucleation, Surface Growth, and Coagulation Problems
A. Prakash , A. P. Bapat & M. R. Zachariah (2003)

Суть: узловой метод решение общего уравнения динамики распределения частиц за счет процессов нуклеации, коагуляции и
гигроскопического роста/уменьшения размеров частиц (eq. 1).

Узловой метод менее затратный, чем метод деления диапазона на непересекающиеся интервалы.
В работе проведено тестирование модели в различных режимах, которое показало хорошее соответствие с имеющимися данными.

Athor: J. Kuznetsova
Date : 07.07.2022
"""

from timeit import default_timer as timer  # загрузка функции-таймера
import os  # модуль для проверки информации о файлах
import Global_parametrs as gp  # глобальные переменные (Входные данные)
import numpy as np  # подключение библиотеки numpy
from math import sqrt, exp
import plotResults
from Coagulation import Coagulation
from Coagulation_Nucleation import Coagulation_Nucleation

start = timer()  # запуск таймера

path_main = "main.inp"
if os.path.exists(path_main) == False and os.path.getsize(path_main) <= 0:
    print("\033[1;35mWarning: \033[1;30m Could not open file!")
    exit()
else:
    fptr = open(path_main, "r+")  # Opening the input file that contains property data for the aerosol material.
    MAX = int(fptr.readline().split()[-1])  # Number of nodes
    T = float(fptr.readline().split()[-1])  # Starting Temperature(Kelvin)
    coolrate = float(fptr.readline().split()[-1])  # Cooling Rate(Kelvin / s): 1000
    P = float(fptr.readline().split()[-1])  # Pressure(Atm): 1
    MW = float(fptr.readline().split()[-1])  # Molecular Weight(Kg / mol): 0.02695
    rho = float(fptr.readline().split()[-1])  # Density(Kg / m3): 2700
    tempLine = fptr.readline().split()[-2:]  # Surface Tension in the form A-BT(dyne/cm):
    A = float(tempLine[0][2:])
    B = float(tempLine[1][2:])
    tempLine = fptr.readline().split()[-2:]  # Saturation Vapor Pressure of material in form Ps=exp(C-D/T)P:
    C = float(tempLine[0][2:])
    D = float(tempLine[1][2:])
    choice = int(fptr.readline().split()[-1])  # Enter your choice -
    # 1)Coagulation,2)nucleation+coagulation,3)nucleation + coagulation + surface growth,4)surface growth:
    beta_option = int(fptr.readline().split()[-1])  # Choice of collision frequency function:
    # 1)Free Molecular Regime,2)Fuchs form for free molecular and transition regime:
    tempLine = fptr.readline().split()[-2:]  # Sutherland's constants for carrier gas in the form A*T^1.5/(B+T):
    A_mu = float(tempLine[0][2:])
    B_mu = float(tempLine[1][2:])
    fptr.close()  # Closing the main property data file.
    q = 10.0**(12.0/(MAX - 1))  # Calculation of geometric spacing factor which depends on the number of nodes.
    v1 = MW/(rho*gp.Na)  # Volume of a monomer unit (a molecule); Na is the Avogadro's Number.
    parameters = {'MAX': MAX, 'P': P, 'T': T, 'rho': rho, 'A_mu': A_mu, 'B_mu': B_mu, 'A': A, 'B': B, 'C': C, 'D': D, 'q': q, 'beta_option': beta_option, 'v1': v1, 'coolrate': coolrate}
    if choice == 1:  # PURE COAGULATION
        v, N = Coagulation(parameters)
    if choice == 2:
        v, N = Coagulation_Nucleation(parameters)

plotResults.plot_distribution(v, N)
plotResults.plot_distribution_NodelToPNC(N)

end = timer()  # остановка таймера
print('время выполнения: %.3e с' % (end - start))
