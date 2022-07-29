# //************************************************************
# //             Nucleation + Coagulation
# //************************************************************
import Global_parametrs as gp  # глобальные переменные (Входные данные)
import numpy as np  # подключение библиотеки numpy
from math import sqrt, exp, log


def Coagulation_Nucleation(parameters):
    MAX = parameters['MAX']
    rho = parameters['rho']
    v1 = parameters['v1']
    q = parameters['q']
    T = parameters['T']
    P = parameters['P']
    A = parameters['A']
    B = parameters['B']
    C = parameters['C']
    D = parameters['D']
    beta_option = parameters['beta_option']
    A_mu = parameters['A_mu']
    B_mu = parameters['B_mu']
    coolrate = parameters['coolrate']

    N = np.zeros(MAX)  # Setting initial number of particles = 0 in all the nodes
    v = np.zeros(MAX)
    v[0] = v1  # Assigning first node to monomer. Note that these are not considered as particles
    dp = np.zeros(MAX)
    dp[0] = (6 * v[0] / gp.pi) ** (1 / 3)  # Setting diametr of particle (? why not = d?)
    m = np.zeros(MAX)
    m[0] = v[0] * rho
    for i in range(1, MAX):
        v[i] = v[0] * q ** i
        dp[i] = (6 * v[i] / gp.pi) ** (1 / 3)
        m[i] = v[i] * rho

    X = np.zeros((MAX, MAX, MAX))  # Calculation of size splitting operators.
    for k in range(1, MAX - 1):
        for i in range(0, MAX):
            for j in range(0, MAX):
                # Conditions in parentheses check if the combined volume of colliding particles is between k and k+1.
                if (v[k] <= (v[i] + v[j])) and ((v[i] + v[j]) < v[k + 1]):
                    X[i][j][k] = (v[k + 1] - v[i] - v[j]) / (v[k + 1] - v[k])
                elif (v[k - 1] <= (v[i] + v[j])) and ((v[i] + v[j]) < v[k]):
                    X[i][j][k] = (v[i] + v[j] - v[k - 1]) / (v[k] - v[k - 1])
                else:
                    X[i][j][k] = 0
    t = 0.0  # Initializing time.
    step = 5e-4  # Timestep for integration
    S = 1.001  # Setting initial saturation ratio,
    # a little larger than 1. Hereafter the saturation ratio is determined
    # by cooling or heating rate of the aerosol. A heating rate would require
    # a negative value

    # Calculating surface area of the monomer, given the volume of monomer v1
    temp1 = 6 * v[0] / gp.pi
    s1 = gp.pi * temp1**(2/3)

    # Calculating mass of the monomer unit
    m1 = rho * v[0]

    # Calculating saturation vapor pressure for the Aerosol.
    Ps = exp(13.07-36373.0/T)*101325.0*P

    # Calculating number concentration of monomers at saturation.
    ns = Ps/(gp.kb * T)

    # Calculating monomer concentration ( multiplying saturation ratio by saturation concentration of monomers).
    N[0] = S * ns

    counter = 0  # Temporary variable for printing after a large number of iterations
    print("\n time \t\t Monomer \t\t Jk \t\t S \t\t Dp_mean \t\t k* \t\t Ntot")
    while T > 300:  # Calculating aerosol properties until the system cools down to 27 C.
        sigma = (A - B*T) * 1e-3  # Surface tension
        Ps = exp(C - D/T) * 101325.0 * P  # Saturation pressure
        ns = Ps / (gp.kb*T)  # Saturation concentration of monomers using ideal gas law
        S = N[0] / ns  # Saturation ratio
        theta = (s1 * sigma)/(gp.kb * T)  # Dimensionless surface tension theta
        a = (2*sigma)/(gp.pi*m1)  # Temporary variable
        b = theta - (4 * theta**3)/(27 * log(S)**2)  # Temporary variable.
        Jk = ns**2 * S * v[0] * a**(1/2) * exp(b)  # Nucleation rate using the classical SCC model.
        c = (3/3) * theta/log(S)  # Temporary variable.
        kstar = c**3  # Calculating critical cluster size that will determine the node at which nucleation occurs.
        dpstar = 4 * sigma*v[0]/(gp.kb * T * log(S))  # Size of particles corresponding to the critical node.
        vstar = gp.pi * dpstar**3 / 6  # Volume of particle corresponding to the critical node.
        Ntot = 0.0  # Initializing total number of particles to zero, as there are no particles prior to nucleation.
        for i in range(MAX):
            Ntot = Ntot + N[i]  # Total number of particles (does not include monomers).
        Vtot = 0.0
        for i in range(MAX):
            Vtot = Vtot + N[i]*v[i]  # Total volume of particles. Note that loop runs
            # from i=2 because i=1 corresponds to monomers,
            # which we do not count as particles.
        Vtotal = Vtot + N[0]*v[0]  # Total volume for mass conservation check (includes monomers)
        Vav = Vtot/Ntot  # Average volume of particles ( number average)
        dpav = (6*Vav/gp.pi)**(1/3)  # Volume based mean diameter of particles printing
        # after every 50 times the loop runs.
        counter = counter + 1
        if counter == 50:
            print("\n%e\t%e\t%e\t%e\t%e\t%e\t%e" % (t, N[0], Jk, S, dpav, kstar, Ntot))
            counter = 0
        # Calculation of collision frequency function beta(i,j)=K(i,j)
        K = np.zeros((MAX, MAX))
        Kmin = 1e-9  # setting an arbitrary value for Kmin to start with.
        if beta_option == 1:
            for i in range(MAX):
                for j in range(MAX):
                    temp1 = 1/v[i] + 1/v[j]
                    temp2 = v[i]**(1/3) + v[j]**(1/3)
                    K[i][j] = (3.0/(4.0*gp.pi))**(0.1666666) * (6.0*gp.kb * T / rho)**(1/2) * temp1**(1/2) * temp2**2
                    if K[i][j] < Kmin:
                        Kmin = K[i][j]  # Calculating the smallest collision frequency
                        # function to decide the characteristic coagulation time.
#

        if beta_option == 2:
            mu = A_mu * T**(3/2)/(B_mu+T)
            lambd = (mu/(P*101325.0)) * sqrt(gp.pi * gp.R * T/(2.0*0.04))
            for i in range(MAX):
                for j in range(MAX):
                    kn1 = (2.0*lambd) / dp[i]
                    kn2 = (2.0*lambd) / dp[j]
                    D1 = (gp.kb * T) / (3.0*gp.pi * mu * dp[i])*((5.0+4.0*kn1+6.0*kn1*kn1+18.0*kn1*kn1*kn1)/(5.0-kn1+(8.0+gp.pi)*kn1*kn1))
                    D2 =(gp.kb * T)/(3.0*gp.pi*mu*dp[j])*((5.0+4.0*kn2+6.0*kn2*kn2+18.0*kn2*kn2*kn2)/(5.0-kn2+(8.0+gp.pi)*kn2*kn2))
                    c1 = sqrt((8.0*gp.kb*T)/(gp.pi*m[i]))
                    c2 = sqrt((8.0*gp.kb*T)/(gp.pi*m[j]))
                    l1 = (8.0*D1)/(gp.pi*c1)
                    l2 = (8.0*D2)/(gp.pi*c2)
                    g1 = ((dp[i]+l1)**3 - (dp[i]*dp[i]+l1*l1)**(3/2)) / (3.0*dp[i]*l1) - dp[i]
                    g2 = ((dp[j]+l2)**3 - (dp[j]*dp[j]+l2*l2)**(3/2)) / (3.0*dp[j]*l2) - dp[j]
                    K[i][j] = 2.0*gp.pi*(D1+D2)*(dp[i]+dp[j])/((dp[i]+dp[j])/(dp[i]+dp[j]+2.0*sqrt(g1*g1+g2*g2))+(8.0*(D1+D2))/(sqrt(c1*c1+c2*c2)*(dp[i]+dp[j])))
                    if K[i][j] < Kmin:
                        Kmin = K[i][j]   # Calculating the smallest collision frequency
                        # function to decide the characteristic coagulation time.

        # Operator to put nucleated particles in the bin just higher than k*
        zeta = np.zeros(MAX)
        for k in range(1, MAX):
            if vstar < v[0]:
                zeta[2] = vstar / v[2]  # Putting particles formed smaller than monomers in the smallest particle node (node 2). This situation arises when k* falls below monomer size.
            elif (v[k-1] <= vstar) and (vstar < v[k]):
                zeta[k] = vstar/v[k]  # Putting particles in node just larger than k*
            else:
                zeta[k] = 0.0

        # Calculation of gain and loss terms for Nk due to coagulation
        for k in range(1, MAX):
            # Initializing gain and loss terms for coagulation
            sum1 = 0.0  # Production term due to collision of two smaller particles to form a k size particle
            sum2 = 0.0  # Loss term due to collision of k size particle with any other particle
            for i in range(1, MAX):
                sum2 = sum2 + K[k][i]*N[i]  # k collides with any other particle to get out of node k, thus loss term
                for j in range(1, k+1):
                    sum1 = sum1 + X[i][j][k] * K[i][j] * N[i]*N[j]  # i and j collide to form particle of size k
                    # which is then multiplied by the size splitting operator to put the particles
                    # in the adjacent nodes after adjusting the volume

            # Change in PSD due to nucleation + coagulation
            N[k] = N[k] + step * (0.5 * sum1 - N[k] * sum2 + Jk * zeta[k])

        # Monomer balance. accounting for the monomers loss due to nucleation.
        N[0] = N[0] - Jk*kstar*step
        T = T - step * coolrate  # Temperature would decrease as time progresses, due to cooling.
        t = t + step  # Time increment

    data = open('res_Coagulation_Nucleation.txt', 'w')
    print("\n\nFinal Particle Size Distribution\n\n")
    print("\nvolume\t\tnumber")
    for i in range(1, MAX):
        print("%e \t %e \n" % (v[i], N[i]))  # Printing the final number distribution obtained due to
        # nucleation and coagulation (number concentration or particles against nodes)
        data.write("%e \t %e \n" % (v[i], N[i]))
    data.close()

    return v, N