import datetime
import matplotlib.pyplot as plt
import numpy as np

# SIR
def ODE_SIR(s0, i0, r0, infR, recR):
    dSdt = -infR * s0 * i0;
    dIdt = infR * s0 * i0 - recR * i0;
    dRdt = recR * i0;
    s1 = s0 + dSdt;
    i1 = i0 + dIdt;
    r1 = r0 + dRdt;
    return s1, i1, r1
    
days = 365

Susceptible = np.zeros(days + 1) 
Infected = np.zeros(days + 1) 
Recovered = np.zeros(days + 1) 


Population = 1000000
InfectionRate = 0.21 
RecoverRate = 0.14 
Infected[0] = 1.0 / Population 
Susceptible[0] = 1.0 - Infected[0] 
Recovered[0] = 0.0 


PeakPatient = 0.0 
PeakDay = 0 
R0 = InfectionRate / RecoverRate 


for t in range(days):
    Susceptible[t + 1], Infected[t + 1], Recovered[t + 1] = ODE_SIR(Susceptible[t], Infected[t], Recovered[t], InfectionRate, RecoverRate)
    if Infected[t + 1] < 1.0 / Population:
        break
    if PeakPatient == 0.0 and Infected[t + 1] < Infected[t]:
        PeakPatient = Infected[t]
        PeakDay = t
    

Infected *= Population
Susceptible *= Population
Recovered *= Population
TerminateDay = t + 1
PeakPatient *= Population
TotalInfection = Recovered[TerminateDay] 


ElapsedDays = (datetime.datetime.today() - datetime.datetime(2020, 1, 16)).days


X = np.zeros(days + 1)
for t in range(days + 1):
    X[t] = t
    
plt.title("SIR Epidemics Model")
plt.xlabel("days")
plt.ylabel("poplation")
plt.plot(X, Susceptible, label="Susceptible")
plt.plot(X, Infected, label = "Infected")
plt.plot(X, Recovered, label = "Recovered")
plt.xlim(0, TerminateDay)


if ElapsedDays < TerminateDay:
    plt.vlines(ElapsedDays, 0, Population, 'red', linestyles='dashed', label='today')
    plt.text(ElapsedDays +2, Population/2, 'Today', rotation=90, verticalalignment='center')

plt.legend()
plt.show()


print("Basic reproduction number (R0): ", R0)
print("Number of infected people: %d/%d (%d%%)" % (TotalInfection, Population, TotalInfection * 100.0 / Population))
print("Peak: %d patients, day %d" % (PeakPatient, PeakDay))
print("Outbreak End: day", TerminateDay)