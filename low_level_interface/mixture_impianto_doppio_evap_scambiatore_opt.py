import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mixture_impianto_doppio_evap_function_T6 import Funz
from scipy.optimize import minimize ,shgo ,dual_annealing,differential_evolution,basinhopping

lib="REFPROP"

fluids="CO2&R1234YF"
#fluids="CO2&PROPANE"
#fluids="CO2&ISOBUTANE"
#fluids="CO2&HEXANE"
#fluids="CO2&R1233ZD"
#fluids="CO2&CO2"






eps=0.8
T_gc=40
P_gc=90
T_eva=-15
T_EVA_AT=0#7
eta_c=0.8
eta1=0.9
eta2=0.8
eta_mix=0.85
v=10 #m/s

T_ref=273.15
dT=10
num=1




mix   = CP.AbstractState(lib, fluids)
mix_l = CP.AbstractState(lib, fluids)
mix_g = CP.AbstractState(lib, fluids)

x=1
#mix.set_mole_fractions([x, 1-x])
mix.set_mass_fractions([x, 1-x])

def f_T_EVA_AT(y):
    f=np.zeros(2)
    
    T_EVA_AT=y[0]
    x=y[1]
    mix.set_mass_fractions([x, 1-x])
    
    funz=Funz(eps,P_gc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
    
    T,P,H,S,m,cop=funz.imp()
    
    xx=(H[7]-H[6])*m[2]/(H[0]-H[9]+(H[7]-H[6])*m[2])
    
    f[0]=(T[9]-(dT*xx+T_eva+T_ref))
    f[1]=(T[7]-(dT*xx+T_eva+T_ref))
    
    print('constr= ',f)
    return f

y=[0,0.95]
y=fsolve(f_T_EVA_AT,y)        
T_EVA_AT=y[0]
x=y[1] 
    
mix.set_mass_fractions([x, 1-x])
funz=Funz(eps,P_gc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
    
T,P,H,S,m,cop=funz.imp()

print(cop)
print(m)



"evaporatore"
n=100
H_e1=np.linspace(H[6],H[7],n)
T_e1=np.zeros(n)
H_e2=np.linspace(H[9],H[0],n)
T_e2=np.zeros(n)
for i in range(n):
    mix.update(CP.HmassP_INPUTS, H_e1[i], P[7])
    T_e1[i]=mix.T()
    mix.update(CP.HmassP_INPUTS, H_e2[i], P[0])
    T_e2[i]=mix.T()

plt.figure(dpi=200)
plt.plot((H_e1-H[6])*m[2]/1000,T_e1-T_ref)
plt.plot((H_e2-H[9]+(H[7]-H[6])*m[2])/1000,T_e2-T_ref)
plt.plot((0,(H[0]-H[9]+(H[7]-H[6])*m[2])/1000),(T_eva,T_eva+dT),'r')
#plt.plot((0,(H[0]-H[9]+(H[7]-H[6])*m[2])/1000),(T[6]+5-T_ref,T[0]+5-T_ref),'r')
plt.xlabel("Q [k]")
plt.ylabel("T [°C]")
plt.title('Evaporatore, cop='+str(np.round(cop,2)))
plt.grid()







