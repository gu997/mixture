import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c
#from mixture_impianto_senza_eiettore_sep_function import Funz
from mixture_impianto_senza_eiettore_sep_function_T7fix import Funz

import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from multiprocessing import Pool
import time


lib="REFPROP"

fluid1="CO2&"
fluid2=[]
fluid2.append("CO2") #cop=2.25
fluid2.append("R1234YF") #cop=2.25
# =============================================================================
# fluid2.append("R1233ZD") #cop=2.74
# fluid2.append("R1234ZEE")#cop=2.43
# fluid2.append("PENTANE" )#cop=2.97
# fluid2.append("HEXANE" ) #cop=3.17
# fluid2.append("PROPANE" ) #cop=3.17
# fluid2.append("ISOBUTANE" ) #cop=3.17
# =============================================================================

fluids_list=[fluid1+e for e in fluid2]


x=0.97
n=20
eps=0.8
T_gc=40
P_gc=np.linspace(70,110,n)#95
cop=np.zeros(n)
T_7=np.zeros(n)
T_8=np.zeros(n)
T_eva=-5
T_sep=10
T_ref=273.15
eta_c=0.8

class F:  
    def __init__(self,mix,mix_l,mix_g):        
        self.mix = mix
        self.mix_l = mix_l
        self.mix_g = mix_g


    def process(self,P_gcc):
        print(P_gcc)
        funz=Funz(eps,*P_gcc,T_gc,T_eva,T_sep,eta_c,self.mix,self.mix_l,self.mix_g)
            
        T,P,H,S,m,cop=funz.imp()
        T,P,H,S=funz.punti_rimanenti(T,P,H,S)
        #print(T[8]-T_ref)
        return cop #, T[7]-T_ref, T[8]-T_ref
  
if __name__ == "__main__":
  
    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True,dpi=200)
    #plt.figure(dpi=200)
    for j,fluids in enumerate(fluids_list):
        print('fluids = '+fluids+'\n'*10)
    
        mix   = CP.AbstractState(lib, fluids)
        mix_l = CP.AbstractState(lib, fluids)
        mix_g = CP.AbstractState(lib, fluids)
        
        f=F(mix,mix_l,mix_g)
        
        #T,P,H,S,m,cop=funz.imp()
        
        #mix.set_mole_fractions([x, 1-x])
        mix.set_mass_fractions([x, 1-x])
        
        params = [[Pi] for Pi in P_gc]
        
        start = time.time()
        cop = list(map(f.process, params)) 
        print(time.time() - start)
        
# =============================================================================
#         start = time.time()
#         with Pool(processes=4) as pool:
#             cop = list(pool.map(f.process, params))    
#         print(time.time() - start)
# =============================================================================

        #cop=Parallel(n_jobs=10)(delayed(process)(i) for i in P_gc)
        
    # =============================================================================
    #     for i in range(n):
    #         print('i =',i)
    #         cop[i],T_7[i],T_8[i]=process(P_gc[i])
    # =============================================================================
    
        
        
        #plt.plot(P_gc,cop,label=fluid2[j])
        ax0.plot(P_gc,cop)
        ax1.plot(P_gc,T_8-T_7,label=fluid2[j])
        #ax1.plot(P_gc,T_8,label=fluid2[j])
    
    
    
    ax1.plot([90,90],[0,0])    
    ax1.set_xlabel("$ P_{gc} $ [bar]")
    ax0.set_ylabel("COP []")
    ax1.set_ylabel("GLIDE [°C]")
    
    ax0.set_title('Curva caratteristica $ P_{gc} $, x='+str(x)+' , $ T_{gc} $=40°C, $ T_{eva} $=-5°C')
    ax0.grid()
    ax1.grid()
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.9),ncol=3)    
    # =============================================================================
    # T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    # 
    # cop2=m[1]*(H[8]-H[7])/(H[2]-H[1])
    # print('cop2=',cop2)
    # 
    # 
    # pc=60*10**5
    # gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g,pc, fluids,x)
    # =============================================================================
