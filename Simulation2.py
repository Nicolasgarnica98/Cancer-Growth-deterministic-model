import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

beta = 0.8 #float(input('Parametro beta: '))
delta = 0.0000864 #float(input('Parametro delta: '))
theta = 0.864 #float(input('Parametro theta: '))
mu = 0.2 #float(input('Parametro mu: '))
gamma = 0.6 #float(input('Parametro gamma: '))
mu_fm = 0.0000005 #float(input('Parametro mu_fm: '))
mig = 0.00002 #float(input('Parametro Mig: '))
K_max = 533000000 #float(input('Máximo numero de CMA: '))
time_max = 1500 #float(input('Insetar cantidad de días para el análisis: '))
def sis_ecuaciones(Ctm_in,t):
    Ctm = Ctm_in[0]
    Ctc = Ctm_in[1]
    migracion = 0
    mult = 1
    if Ctm/K_max <1:
        mult = 0
    migracion = (mig*(Ctm**2)*theta*mu_fm*beta)*mult #*(np.floor(Ctm/K_max))
    dCtmdt = (((mu*Ctm))-(delta*Ctm)-migracion)*(1-int(Ctm/K_max))
    dCtcdt = (((mu-delta)*Ctc) + (gamma*migracion))*(1-int(Ctc/K_max))
    return [dCtmdt,dCtcdt]

Ctm_0 = 1 #float(input('Células inciales en tumor original: ')) #Condiciones iniciales 1
Ctc_0 = 0 #float(input('Células cancerigenas iniciales en tejido objetivo: ')) #Condiciones iniciales 2
t = np.arange(0,int(time_max),0.001) #Variable independiente

sol = odeint(sis_ecuaciones,[Ctm_0,Ctc_0],t)
print(sol)

#Grafica
plt.title('Simulación de proliferación celular')
plt.plot(t,sol[:,0],label= "Tumor pulmonar Metástico: Ctm(t)", color='green')
plt.plot(t,sol[:,1],label= "Tumor cerebral objetivo: Ctc(t)", color='red')
plt.hlines(K_max,0,time_max,linestyles="--", label= "No. máximo de CACs")
plt.legend()
plt.ylabel("No. Células")
plt.xlabel("Tiempo (Días)")
plt.tight_layout()
plt.grid()
plt.show()