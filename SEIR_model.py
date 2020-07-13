import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 10000000

# Initial number of infected and recovered individuals, I0 and R0.
E0,I0, R0 = 249*2.399, 563, 53
# Everyone else, S0, is susceptible to infection initially.
S0 = N - E0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
for iteration in range(0,11,1):
    rho=iteration/10
    beta, gamma, sigma = rho*2.75/2.3, 1./2.3, 1./5.2  
    # A grid of time points (in days)
    t = np.linspace(0, 400, 400)

    # The SEIR model differential equations.
    def deriv(y, t, N, beta, gamma,sigma):
        S, E, I, R = y
        dSdt = -beta * S * I / N
        dEdt=  beta*S*I/N - sigma*E
        dIdt = (sigma * E- gamma* I)
        dRdt = gamma * I
        return dSdt, dEdt, dIdt, dRdt

    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # Integrate the SEIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, sigma))
    S, E, I, R = ret.T

    # Plot the data on three separate curves for S(t), E(t), I(t) and R(t)
    fig = plt.figure(facecolor='w',figsize=(8,8))
    ax = fig.add_subplot(111, facecolor='w',axisbelow=True)
    #ax.plot(t, S/1000000, '#3792cb', alpha=1, lw=2, label='Susceptible')
    ax.plot(t, E/100000000, '#ffaa1d', alpha=1, lw=2, label='Exposed')
    ax.plot(t, I/100000000, '#d21f3c', alpha=1, lw=2, label='Infected')
    #ax.plot(t, R/1000000, '#0b6623', alpha=1, lw=2, label='Recovered')
    ax.set_xlabel('Days')
    ax.set_ylabel('Population (10s of crores)')
    ax.set_ylim(0,2.5)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    plt.title("COVID-19 SEIR model showcasing the impact of social distancing")

    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)

    textstr = '\n'.join((
        r'$\rho=%.2f$' % (rho, ),
        r'$\beta=%.2f$' % (beta, ),
        r'$\gamma=%.2f$' % (gamma, ),
        r'$\sigma=%.2f$' % (sigma, )))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='w', alpha=1)

    # place a text box in upper left in axes coords
    ax.text(0.92, 0.85, textstr, transform=plt.gcf().transFigure, fontsize=10,
            verticalalignment='top', bbox=props)
    fig1 = plt.gcf()
    plt.show()