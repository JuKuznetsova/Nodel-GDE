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




    # //**********************************************
    # //           PURE COAGULATION
    # //**********************************************
    # This part of the code solves for coagulation problems.
    # The code verifies that we obtain a self preserving distribution.
    # It verifies the problem stated in Friedlander, 2000 - page 218, line 4 - and shows that SPD is obtained.
    # And note that in this case node 1 is not the monomer. Its particle of size 1 nm.
    # We do not have to worry about monomers, since we are concerned with coagulation only.

    if choice == 1:
        path_coag = "coag.inp"
        if os.path.exists(path_coag)==False and os.path.getsize(path_coag) <= 0:
            print("\033[1;35mWarning: \033[1;30m Could not open file!")
            exit()
        else:
            fptr1 = open(path_coag, "r+")  # Opening the input file that contains property data for the aerosol material.
            N = np.zeros(MAX)

            # Scanning the initial monomer size and number concentration of the monodisperse distribution for coagulation.
            d = float(fptr1.readline().split()[-1])  # Initial diameter of monodisperse particles(nm):
            N[0] = float(fptr1.readline().split()[-1])  # Initial number concentration(#/m3):
            T = float(fptr1.readline().split()[-1])  # Temperature(Kelvin):
            fptr1.close()  # Closing the coag property data

            v = np.zeros(MAX)
            v[0] = 1e-27 * gp.pi * d**3 / 6  # Setting volume of particle
            # that the user enters as the smallest node(coagulation would never lead to decrease in particle size).
            dp = np.zeros(MAX)
            dp[0] = (6*v[0]/gp.pi)**(1/3)  # Setting diametr of particle (? why not = d?)
            m = np.zeros(MAX)
            m[0] = v[0] * rho
            for i in range(1, MAX):
                v[i] = v[0] * q**i
                dp[i] = (6 * v[i] / gp.pi)**(1/3)
                m[i] = v[i] * rho

            X = np.zeros((MAX, MAX, MAX))  # Calculation of size splitting operators.
            for k in range(1, MAX-1):
                for i in range(0, MAX):
                    for j in range(0, MAX):
                        # Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
                        if (v[k] <= (v[i] + v[j]) and (v[i] + v[j]) < v[k + 1]):
                            X[i][j][k] = (v[k + 1] - v[i] - v[j]) / (v[k + 1] - v[k])
                        elif (v[k-1] <= (v[i] + v[j]) and (v[i] + v[j]) < v[k]):
                            X[i][j][k] = (v[i] + v[j] - v[k - 1]) / (v[k] - v[k - 1])
                        else:
                            X[i][j][k] = 0

            t = 0.0  # Initializing time.
            step = 1e-12  # Initializing timestep for integration.
            t_SPD = 5 / ((3.0 / (4.0 * gp.pi))**(1/6) * (6.0 * gp.kb * T / rho)**0.5 * v[0]**(1/6) * N[0])
            print("\n\nEstimated time required to reach SPD = %e" % t_SPD)
            Ninf_eq = 1e24  # A temporary variable that calculates N_total according to Friedlander, 2000, equation 7.77.

            print("\n t \t\t Ninf \t\t N[1] \t\t N[5] \t\t N[10] \t\t N[15]")
            # Running the coagulation mode of the code until SPD is reached.
            # The dimensionless number distribution remains same after SPD has reached.
            K = np.zeros((MAX, MAX))
            while (t < 1e-2):  #(t < t_SPD):
                # Calculation of collision frequency function ( same as in " nucleation + coagulation" section of the code).
                Kmin = 1e-9  # setting an arbitrary value for Kmin to start with.
                if beta_option == 1:
                    for i in range(0, MAX):
                        for j in range(0, MAX):
                            temp1 = 1 / v[i] + 1 / v[j]
                            temp2 = v[i]**(1/3) + v[j]**(1/3)
                            K[i][j] = (3.0 / (4.0 * gp.pi))**(1/6) * (6.0 * gp.kb * T / rho)**(1/2) * temp1**(1/2) * temp2**2
                            if (K[i][j] < Kmin):  # Calculating the smallest collision frequency function to decide the characteristic coagulation time.
                                Kmin = K[i][j]
                if beta_option == 2:
                    mu = A_mu * T**1.5 / (B_mu + T)
                    lambd = (mu / (P * 101325.0)) * sqrt(gp.pi * gp.R * T / (2.0 * 0.04))
                    for i in range(0, MAX):
                        kn1 = (2.0 * lambd) / dp[i]
                        for j in range(0, MAX):
                            kn2 = (2.0 * lambd) / dp[j]
                            D1 = (gp.kb * T) / (3.0 * gp.pi * mu * dp[i]) * (
                                        (5.0 + 4.0 * kn1 + 6.0 * kn1 * kn1 + 18.0 * kn1 * kn1 * kn1) / (
                                            5.0 - kn1 + (8.0 + gp.pi) * kn1 * kn1))
                            D2 = (gp.kb * T) / (3.0 * gp.pi * mu * dp[j]) * (
                                        (5.0 + 4.0 * kn2 + 6.0 * kn2 * kn2 + 18.0 * kn2 * kn2 * kn2) / (
                                            5.0 - kn2 + (8.0 + gp.pi) * kn2 * kn2))
                            c1 = sqrt((8.0 * gp.kb * T) / (gp.pi * m[i]))
                            c2 = sqrt((8.0 * gp.kb * T) / (gp.pi * m[j]))
                            l1 = (8.0 * D1) / (gp.pi * c1)
                            l2 = (8.0 * D2) / (gp.pi * c2)
                            g1 = ((dp[i] + l1)**3 - (dp[i] * dp[i] + l1 * l1)**1.5) / (3.0 * dp[i] * l1) - dp[i]
                            g2 = ((dp[j] + l2)**3 - (dp[j] * dp[j] + l2 * l2)**1.5) / (3.0 * dp[j] * l2) - dp[j]
                            K[i][j] = 2.0 * gp.pi * (D1 + D2) * (dp[i] + dp[j]) / (
                                        (dp[i] + dp[j]) / (dp[i] + dp[j] + 2.0 * sqrt(g1 * g1 + g2 * g2)) + (
                                            8.0 * (D1 + D2)) / (sqrt(c1 * c1 + c2 * c2) * (dp[i] + dp[j])))
                            # print("kn %e\n", kn2)
                            if (K[i][j] < Kmin):  # Calculating the smallest collision frequency function to decide the characteristic coagulation time.
                                Kmin = K[i][j]

                for k in range(0, MAX):  # Calculating the gain and loss terms due to coagulation.
                    sum1 = 0  # Addition term when i and j collide to form a k sized particle.
                    sum2 = 0  # Subtraction term when k collides with any other particle.
                    for i in range(0, MAX):
                        sum2 = sum2 + K[k][i] * N[i]
                        for j in range(0, k+1):
                            sum1 = sum1 + X[i][j][k] * K[i][j] * N[i] * N[j]
                    N[k] = N[k] + step * (0.5 * sum1 - N[k] * sum2)

                Ninf = 0  # Total number of particles, which is zero initially.
                for i in range(0, MAX):
                    Ninf = Ninf + N[i]  # N_infinity as we get by solving the coagulation part of GDE.

                step = 1e-3 / (Kmin * Ninf)  # Adaptive timestep for integration, based on characteristic coagulation time.
                fi = 0
                for i in range(1, MAX):
                    fi = fi + N[i] * v[i]
                # N_infinity according to equation 7.77 in Friedlander
                # Ninf_eq = Ninf_eq - step * (3.33 * (3 / (4 * gp.pi))**(1.6) * (6 * gp.kb * T / rho)**(1.2) * fi**(1.6) * Ninf_eq**(11.0 / 6.0))
                Vtot = 0
                print("\n %e \t %e \t %e \t %e \t %e \t %e" % (t, Ninf, N[1], N[5], N[10], N[15]))
                t = t + step

            data = open('res_Coagulation_0.01s.txt', 'w')
            # Printing the final number distribution after SPD is reached.
            print("\n\n***Size Distribution after reaching SPD****\n")
            print("\n volume \t\t number")
            for i in range(1, MAX):
                print("%e \t %e" % (v[i], N[i]))
                data.write("%e \t %e \n" % (v[i], N[i]))
            data.close()

plotResults.plot_distribution(v, N)


end = timer()  # остановка таймера
print('время выполнения: %.3e с' % (end - start))