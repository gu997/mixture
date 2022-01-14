import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c
import matplotlib.pyplot as plt


class Funz:  
    def __init__(self,eps,P_gc,T_gc,T_eva,T_sep,eta_c,mix,mix_l,mix_g):        
        

        self.mix = mix
        self.mix_l = mix_l
        self.mix_g = mix_g

        self.eps=eps
        self.T_gc=T_gc
        self.P_gc=P_gc
        self.T_eva=T_eva
        self.T_sep=T_sep
        self.eta_c=eta_c
       
        self.T=np.zeros(11)
        self.P=np.zeros(11)
        self.H=np.zeros(11)
        self.S=np.zeros(11)
        self.m=np.ones(3)
        
        
        self.T_ref=273.15
               


    def imp(self):
        "PUNTI NOTI A PRIORI"
        
        "punto 3"
        self.T[3]=self.T_gc+self.T_ref
        self.P[3]=self.P_gc*10**5   #iperparametro
        self.mix.update(CP.PT_INPUTS, self.P[3], self.T[3])
        self.H[3]=self.mix.hmass()
        self.S[3]=self.mix.smass()
        
        "punto 8"
        self.T[8]=self.T_eva+self.T_ref
                
        "punto 5"
        self.T[5]=self.T_sep+self.T_ref
        
        "punto 4"
        #self.P[4]=self.P[3]
        
        "punto 2"
        self.P[2]=self.P[3]
        
        "IPOTIZZO IL titolo in  5"
        def loop(Q):
            if Q<0:
                Q=0
            
            "punto 5"
            #Q=0.3 #molare
            self.mix.update(CP.QT_INPUTS, Q, self.T[5])
            self.P[5]=self.mix.p()
            self.H[5]=self.mix.hmass()
            
            "SEPARATORE"
            m_l=self.mix.mole_fractions_liquid()
            m_g=self.mix.mole_fractions_vapor()
            
            
            self.mix_l.set_mole_fractions(m_l)
            self.mix_g.set_mole_fractions(m_g)
            
            "PUNTO 8"
            self.mix_l.update(CP.QT_INPUTS, 1, self.T[8])
            self.P[8]=self.mix_l.p()
            self.H[8]=self.mix_l.hmass()
            
            "PUNTO 9"
            self.T[9]=self.T[5]
            self.P[9]=self.P[5]
            self.mix_g.update(CP.PT_INPUTS, self.P[9], self.T[9])
            #mix_g.Q()
            self.H[9]=self.mix_g.hmass()
            
            "PUNTO 10"
            self.H[10]=self.H[9]
            #self.P[10]=self.P[8]
            #mix_l.update(CP.HmassP_INPUTS, H[10], P[10])
            #T[10]=mix_g.T()
            
            "PUNTO 0"
            q=Q*self.mix_g.molar_mass()/self.mix.molar_mass() #massico
            self.m[1]=1-q
            self.m[2]=q
            self.H[0]=self.m[1]*self.H[8]+self.m[2]*self.H[10]
            self.P[0]=self.P[8]
            self.mix.update(CP.HmassP_INPUTS, self.H[0], self.P[0])
            self.T[0]=self.mix.T()
            #S[0]=mix.smass()
            
            "punto 1"
            self.P[1]=self.P[8]
            self.T[1]=self.T[0] + self.eps*(self.T[3]-self.T[0])
            self.mix.update(CP.PT_INPUTS, self.P[1], self.T[1])
            self.H[1]=self.mix.hmass()
            
            "PUNTO 4"
            self.H[4]=self.H[3]-self.H[1]+self.H[0]
            
            "punto 5"
            self.H[5]=self.H[4]
            self.mix.update(CP.HmassT_INPUTS, self.H[5], self.T[5])
        
           
            return self.mix.Q()-Q
        
        Q=fsolve(loop,0.3)
  
        
        "punto 1"
        self.mix.update(CP.PT_INPUTS, self.P[1], self.T[1])
        self.S[1]=self.mix.smass()
        
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
    
        return self.T,self.P,self.H,self.S,self.m,cop



    def punti_rimanenti(self,T,P,H,S):    
        
        
        "PUNTO 6"
        P[6]=P[5]
        T[6]=T[5]
        self.mix_l.update(CP.PT_INPUTS, P[6], T[6])
        #mix_l.Q()
        H[6]=self.mix_l.hmass()
        S[6]=self.mix_l.smass()
        
        
        "PUNTO 7"
        H[7]=H[6]
        P[7]=P[8]
        self.mix_l.update(CP.HmassP_INPUTS, H[7], P[7])
        T[7]=self.mix_l.T()
        S[7]=self.mix_l.smass()
        
        
        
        return T,P,H,S
    
    def entropie(self,T,P,H,S):    
        "PUNTO 0"
        self.mix.update(CP.HmassP_INPUTS, self.H[0], self.P[0])
        S[0]=self.mix.smass()
        
        "punto 5"
        self.mix.update(CP.HmassT_INPUTS, self.H[5], self.T[5])
        S[5]=self.mix.smass()
        
        "PUNTO 4"
        P[4]=P[3]
        self.mix.update(CP.HmassP_INPUTS, H[4], P[4])
        T[4]=self.mix.T()
        S[4]=self.mix.smass()
        
        "PUNTO 8"
        self.mix_l.update(CP.QT_INPUTS, 1, self.T[8])
        S[8]=self.mix_l.smass()
        
        "PUNTO 9"
        self.mix_g.update(CP.PT_INPUTS, self.P[9], self.T[9])
        S[9]=self.mix_g.smass()
        
        
        "PUNTO 10"
        P[10]=P[8]
        self.mix_g.update(CP.HmassP_INPUTS, H[10], P[10])
        T[10]=self.mix_g.T()
        S[10]=self.mix_g.smass()
        
        return T,P,H,S

if __name__ == "__main__":
    
    lib="REFPROP"

    #fluids="CO2&R1234YF"
    #fluids="CO2&PROPANE"
    #fluids="CO2&ISOBUTANE"
    #fluids="CO2&HEXANE"
    #fluids="CO2&R1233ZD"
    #fluids="CO2&methanol" #dà troppo glide e Pcrit elevata
    #fluids="CO2&ethanol" #curva fatta male converge con 0.998
    #fluids="CO2&dimethylether" #sembra buono con x=0.97 
    fluids="CO2&R32"




    mix   = CP.AbstractState(lib, fluids)
    mix_l = CP.AbstractState(lib, fluids)
    mix_g = CP.AbstractState(lib, fluids)
    
    x= 0.7#0.8904072#0.890685902909111#0.8906859#0.94470172#0.99
    #mix.set_mole_fractions([x, 1-x])
    mix.set_mass_fractions([x, 1-x])
    print('P_crit=',mix.keyed_output(CP.iP_critical)/100000)
    
    eps=0.8
    T_gc=40
    P_gc=90#74.21052631578948#88.94736842105263#70#95
    T_eva=-5
    T_sep=10
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
    funz=Funz(eps,P_gc,T_gc,T_eva,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    
    print(cop)
    
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    T,P,H,S=funz.entropie(T,P,H,S)

    pc=60*10**5
    gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g,pc, fluids,x)
    #gt.grafico_TS_sep_IHX(T,S/1000, mix, mix_l, mix_g,pc, fluids,x)
    
    "evaporatore"
    n=100
    H_e=np.linspace(H[7],H[8],n)
    T_e=np.zeros(n)
    for i in range(n):
        mix_l.update(CP.HmassP_INPUTS, H_e[i], P[7])
        T_e[i]=mix_l.T()
    
    plt.figure(dpi=200)
    plt.plot(H_e/1000,T_e-T_ref)
    plt.plot((H[7]/1000,H[8]/1000),(T[7]+5-T_ref,T[8]+5-T_ref),'r')
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("T [°C]")
    plt.title('Evaporatore')
    plt.grid()