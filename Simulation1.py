
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

beta = float(input('Parametro beta: '))
delta = float(input('Parametro delta: '))
theta = float(input('Parametro theta: '))
mu = float(input('Parametro mu: '))
gamma = float(input('Parametro gamma: '))

def sis_ecuaciones(Ctm_in,t):
    Ctm = Ctm_in[0]
    Ctc = Ctm_in[1]
    dCtmdt = (mu-delta-(beta*theta))*Ctm
    dCtcdt = ((mu-delta)*Ctc)+(beta*gamma*theta*Ctm)
    return [dCtmdt,dCtcdt]

Ctm_0 = float(input('Células inciales en tumor original: ')) #Condiciones iniciales 1
Ctc_0 = float(input('Células cancerigenas iniciales en tejido objetivo: ')) #Condiciones iniciales 2
t = np.arange(0,3,0.001) #Variable independiente

sol = odeint(sis_ecuaciones,[Ctm_0,Ctc_0],t)
print(sol)

#Grafica
plt.title('Simulación de proliferación celular')
plt.plot(t,sol[:,0],label= "Tumor original: Ctm(t)", color='green')
plt.plot(t,sol[:,1],label= "Tejido objetivo: Ctc(t)", color='red')
plt.legend()
plt.ylabel("No. Células")
plt.xlabel("Tiempo (Días)")
plt.tight_layout()
plt.grid()
plt.show()