import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mixture_impianto_doppio_evap_function import Funz
from scipy.optimize import minimize ,shgo ,dual_annealing,differential_evolution,basinhopping

lib="REFPROP"

fluids="CO2&R1234YF"
#fluids="CO2&PROPANE"
#fluids="CO2&ISOBUTANE"
#fluids="CO2&HEXANE"
#fluids="CO2&R1233ZD"
#fluids="CO2&CO2"



mix   = CP.AbstractState(lib, fluids)
mix_l = CP.AbstractState(lib, fluids)
mix_g = CP.AbstractState(lib, fluids)

x= 0.85
#mix.set_mole_fractions([x, 1-x])
mix.set_mass_fractions([x, 1-x])


eps=0.8
T_gc=40
P_gc=90
T_eva=-10
T_EVA_AT=0#7
eta_c=0.8
eta1=0.9
eta2=0.8
eta_mix=0.85
v=10 #m/s

T_ref=273.15
dT=10
num=1

funz=Funz(eps,P_gc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
    
T,P,H,S,m,cop=funz.imp()

print(cop)


#T,P,H,S=funz.punti_rimanenti(T,P,H,S)
#T,P,H,S=funz.entropie(T,P,H,S)

cop2=(m[2]*(H[7]-H[6])+m[0]*(H[0]-H[9]))/(H[2]-H[1])/m[0]
print(cop2)
#print(T[8]-T[7])
# =============================================================================
#     pc=60*10**5
#     gt.grafico_PH_eiettore_sep(P/100000,H/1000, mix, mix_l, mix_g,pc, "CO2&R1234YF",x)
# =============================================================================
class arg:  
    def __init__(self,H,T,m): 
        self.H=H
        self.T=T
        self.m=m
      
    def f_glide(self,y):
        print(y)
        T_eva=y[0]*num
        T_EVA_AT=y[1]*num
        x= y[2]
        try:
            mix.set_mass_fractions([x, 1-x])
            funz=Funz(eps,P_gc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
        
            self.T,P,self.H,S,self.m,cop=funz.imp()
            #print(self.T)
        except Exception as e:
            print('\n'*10)
            print(e)
            print('\n'*10)
            cop=0
        return -cop
    def con1(self,y):
        print('con1= ',-(y[1]-T_EVA_AT))
        return -(y[1]*num-T_EVA_AT)
    def con2(self,y):
        print('con2= ',-(self.T[6]-T_EVA_AT+dT-T_ref))
        return -(self.T[6]-T_EVA_AT+dT-T_ref)
    def con3(self,y):
        xx=(self.H[7]-self.H[6])*self.m[2]/(self.H[0]-self.H[9]+(self.H[7]-self.H[6])*self.m[2])
        print('con3= ',-(self.T[7]-(dT*xx+T_EVA_AT-dT+T_ref)))
        return -(self.T[7]-(dT*xx+T_EVA_AT-dT+T_ref))
    def con4(self,y):
        xx=(self.H[7]-self.H[6])*self.m[2]/(self.H[0]-self.H[9]+(self.H[7]-self.H[6])*self.m[2])
        print('con4= ',-(self.T[9]-(dT*xx+T_EVA_AT-dT+T_ref)))
        return -(self.T[9]-(dT*xx+T_EVA_AT-dT+T_ref))

H=np.zeros(10)
T=np.zeros(10)
m=np.zeros(3)
args=arg(H,T,m)
cons = ({'type': 'ineq', 'fun': args.con1},
        {'type': 'ineq', 'fun': args.con2},
        {'type': 'ineq', 'fun': args.con3},
        {'type': 'ineq', 'fun': args.con4})

bnds = ((-20/num, -5/num), (-10/num, 0),(0.7,0.99))
y=[(T_eva+4)/num,(T_EVA_AT-2)/num,x]


res=minimize(args.f_glide,y,method='SLSQP', bounds=bnds, constraints=cons, tol=1e-2)    
#res = shgo(args.f_glide, bounds=bnds,n=64, constraints=cons)


y=res.x



#mix.set_mole_fractions([x, 1-x])
mix.set_mass_fractions([y[2], 1-y[2]])

funz=Funz(eps,P_gc,T_gc,y[0]*num,y[1]*num,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
    
T,P,H,S,m,cop=funz.imp()


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
plt.plot((0,(H[0]-H[9]+(H[7]-H[6])*m[2])/1000),(T_EVA_AT-dT,T_EVA_AT),'r')
#plt.plot((0,(H[0]-H[9]+(H[7]-H[6])*m[2])/1000),(T[6]+5-T_ref,T[0]+5-T_ref),'r')
plt.xlabel("H [kJ/kg]")
plt.ylabel("T [°C]")
plt.title('Evaporatore')
plt.grid()







