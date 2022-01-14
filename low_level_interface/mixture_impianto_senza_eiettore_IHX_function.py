import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici as gt
from scipy.optimize import fsolve
import compressore as c
import matplotlib.pyplot as plt


class Funz:  
    def __init__(self,eps,P_gc,T_gc,T_eva,eta_c,mix,mix_l,mix_g):        
        

        self.mix = mix
        self.mix_l = mix_l
        self.mix_g = mix_g

        self.eps=eps
        self.T_gc=T_gc
        self.P_gc=P_gc
        self.T_eva=T_eva
        self.eta_c=eta_c
       
        self.T=np.zeros(6)
        self.P=np.zeros(6)
        self.H=np.zeros(6)
        self.S=np.zeros(6)
        
        
        self.T_ref=273.15
               


    def imp(self):
        "PUNTI NOTI A PRIORI"
        
        "punto 3"
        self.T[3]=self.T_gc+self.T_ref
        self.P[3]=self.P_gc*10**5   #iperparametro
        self.mix.update(CP.PT_INPUTS, self.P[3], self.T[3])
        self.H[3]=self.mix.hmass()
        self.S[3]=self.mix.smass()
        
        "punto 0"
        self.T[0]=self.T_eva+self.T_ref
        self.mix.update(CP.QT_INPUTS, 1, self.T[0])
        self.P[0]=self.mix.p()
        self.H[0]=self.mix.hmass()
        
        "punto 2"
        self.P[2]=self.P[3]
       
       
        
        "punto 1"
        self.P[1]=self.P[0]
        self.T[1]=self.T[0] + self.eps*(self.T[3]-self.T[0])
        self.mix.update(CP.PT_INPUTS, self.P[1], self.T[1])
        self.H[1]=self.mix.hmass()
        self.S[1]=self.mix.smass()
        
        "PUNTO 4"
        self.P[4]=self.P[3]
        self.H[4]=self.H[3]-self.H[1]+self.H[0]
        self.mix.update(CP.HmassP_INPUTS, self.H[4], self.P[4])
        self.T[4]=self.mix.T()
        
        "punto 5"
        self.P[5]=self.P[0]
        self.H[5]=self.H[4]
        self.mix.update(CP.HmassP_INPUTS, self.H[5], self.P[5])
        self.T[5]=self.mix.T()
        
        "PUNTO 2"
        self.mix.update(CP.PSmass_INPUTS, self.P[2], self.S[1])
        H2_id=self.mix.hmass()
        self.H[2]=self.H[1]+(H2_id-self.H[1])/self.eta_c
        self.mix.update(CP.HmassP_INPUTS, self.H[2], self.P[2])
        self.T[2]=self.mix.T()
        self.S[2]=self.mix.smass()
  
        
# =============================================================================
#         "punto 2"
#         self.T[2]=c.Temperatura_mandata(self.T[0]-self.T_ref, self.T[1]-self.T_ref, self.P[2]*10**-5)+self.T_ref
#         self.mix.update(CP.PT_INPUTS, self.P[2], self.T[2])
#         self.H[2]=self.mix.hmass()
# =============================================================================
        
        
        #cop2=self.m[1]*(self.H[8]-self.H[7])/(self.H[2]-self.H[1])
        cop=(self.H[1]-self.H[3])/(self.H[2]-self.H[1])
        #print('cop =',cop)
        #print('cop2=',cop2)
    
        return self.T,self.P,self.H,self.S,cop



 
if __name__ == "__main__":
    
    lib="REFPROP"

    fluids="CO2&R1234YF"
    #fluids="CO2&PROPANE"
    #fluids="CO2&ISOBUTANE"
    #fluids="CO2&HEXANE"
    #fluids="CO2&R1233ZD"



    mix   = CP.AbstractState(lib, fluids)
    mix_l = CP.AbstractState(lib, fluids)
    mix_g = CP.AbstractState(lib, fluids)
    
    x= 0.8#0.8904072#0.890685902909111#0.8906859#0.94470172#0.99
    #mix.set_mole_fractions([x, 1-x])
    mix.set_mass_fractions([x, 1-x])
    
    
    eps=0.8
    T_gc=40
    P_gc=90#74.21052631578948#88.94736842105263#70#95
    T_eva=-5
    eta_c=0.9
    
    
# =============================================================================
#     T=np.zeros(11)
#     P=np.zeros(11)
#     H=np.zeros(11)
#     S=np.zeros(11)
#     m=np.ones(3)
#     m_ll=np.zeros(2)
#     m_gg=np.zeros(2)
# =============================================================================
    
    
    T_ref=273.15
    
# =============================================================================
#     funz=Funz(eps,P_gc,T_gc,T_eva,T_sep,mix,mix_l,mix_g)
# 
#     funz.imp()
# =============================================================================
    funz=Funz(eps,P_gc,T_gc,T_eva,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,cop=funz.imp()
    
    print(cop)
    
   

    #pc=60*10**5
    plt.figure(dpi=200)
    gt.grafico_PH_semplice(P/100000,H/1000, 'r',1)
    #gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g,pc, fluids,x)
    #gt.grafico_TS_sep_IHX(T,S/1000, mix, mix_l, mix_g,pc, fluids,x)
    
    "evaporatore"
    n=100
    H_e=np.linspace(H[5],H[0],n)
    T_e=np.zeros(n)
    for i in range(n):
        mix.update(CP.HmassP_INPUTS, H_e[i], P[5])
        T_e[i]=mix.T()
    
    plt.figure(dpi=200)
    plt.plot(H_e/1000,T_e-T_ref)
    plt.plot((H[5]/1000,H[0]/1000),(T[5]+5-T_ref,T[0]+5-T_ref),'r')
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("T [°C]")
    plt.title('Evaporatore')
    plt.grid()