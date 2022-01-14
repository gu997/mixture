import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
#import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
from mixture_impianto_senza_eiettore_sep_function_T7fix import Funz as Funz7
from mixture_impianto_senza_eiettore_sep_function import Funz
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#from joblib import Parallel, delayed



lib="REFPROP"

fluid1="CO2&"
fluid2=[] ; x0=[]

fluid2.append("R1234YF") ; x0.append(0.9112132132132) 
# =============================================================================
# fluid2.append("R1233ZD") ; x0.append(0.98) 
# #fluid2.append("R1233ZD") ; x0.append(0.985) #se T7=-20
# 
# 
# fluid2.append("R1234ZEE"); x0.append(0.98) 
# fluid2.append("PENTANE" ); x0.append(0.9955) 
# #fluid2.append("HEXANE" ) ; x0.append(0.99) 
# fluid2.append("PROPANE" ) ; x0.append(0.9) 
# fluid2.append("ISOBUTANE" ) ; x0.append(0.97) 
# fluid2.append("CO2"); x0.append(0.99)  #cop=2.25
# =============================================================================

fluids_list=[fluid1+e for e in fluid2]
fluids=fluids_list[0]


n=200
eps=0.8
T_gc=40
#10**np.linspace(0,np.log10(100),10)    =   np.logspace(0,np.log10(100),10)
P_gc=90#np.linspace(70,110,n)#95
#P_gc=np.linspace(110,70,n)#95
cop=np.zeros(n)
x=np.linspace(0.85,0.9,n)
T_eva=-5
T_sep=10
T_ref=273.15
eta_c=0.8
xx=[]
copcop=[]





glide=10

def process(P_gcc):
    funz=Funz(eps,P_gcc,T_gc,T_eva,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    print(T[8]-T_ref)
    return cop ,  T[8]-T[7]

def process7(P_gcc):
    funz=Funz7(eps,P_gcc,T_gc,T_eva-glide,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    print(T[8]-T_ref)
    return cop ,  T[8]-T[7]
    
def f_glide(x,args):

    P_gcc=args
    mix.set_mass_fractions([x, 1-x])
    cop,dT=process(P_gcc)
    print(x,dT)
    if dT<glide:
        cop,dT=process7(P_gcc)
        print(dT)
    return cop#,dT-glide
    
    

mix   = CP.AbstractState(lib, fluids)
mix_l = CP.AbstractState(lib, fluids)
mix_g = CP.AbstractState(lib, fluids)
    
for i in range(n):
    print('i =',i)
    cop[i]=f_glide(x[i],P_gc)
    
  
plt.figure(dpi=200)
plt.plot(x,cop)
plt.grid()
   
