import scipy.integrate
import numpy
import matplotlib.pyplot as plt
import math 

def SIR_model_WD(y, t, beta, gamma, K, mu):
    S, I , D,  R = y
    
    dS_dt = -K*S*math.log(1+ (beta*I/K))
    dI_dt = K*S*math.log(1+ (beta*I/K)) - (gamma + mu)*I 
    dR_dt = gamma*I
    dD_dt = mu*I
    
    return ([dS_dt, dI_dt, dR_dt, dD_dt])
# Initial conditions

S0 = 9999
I0 = 1
R0 = 0.0 
D0 = 0.0 
beta = 0.9
gamma = 0.15
K = 0.010
mu = 0.05

#time vector

t = numpy.linspace(0, 40, 10000)

#solution

solution = scipy.integrate.odeint(SIR_model_WD, [S0, I0, R0, D0], t , args = (beta, gamma, K, mu))
solution = numpy.array(solution)
#Result

plt.figure(figsize = [8,4])
plt.plot(t, solution[:, 0], label = "S(t)")
plt.plot(t, solution[:, 1], label = "I(t)")
plt.plot(t, solution[:, 2], label = "R(t)")
plt.plot(t, solution[:, 3], label = "D(t)")
plt.xlabel("Time")
plt.ylabel("Individuals (in %)")
plt.grid()
plt.show() 
