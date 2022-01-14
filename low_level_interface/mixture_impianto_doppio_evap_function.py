import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


class Funz:  
    def __init__(self,eps,P_gc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g):        
        

        self.mix = mix
        self.mix_l = mix_l
        self.mix_g = mix_g

        self.eps=eps
        self.T_gc=T_gc
        self.P_gc=P_gc
        self.T_eva=T_eva
        self.T_EVA_AT=T_EVA_AT
        self.eta_c=eta_c
        self.eta1=eta1
        self.eta2=eta2
        self.eta_mix=eta_mix
        self.v=v
       
        self.T=np.zeros(10)
        self.P=np.zeros(10)
        self.H=np.zeros(10)
        self.S=np.zeros(10)
        self.m=np.ones(3)
        
        
        self.T_ref=273.15
               


    def imp(self):
        
        "PUNTI NOTI A PRIORI" #non ci sono
        "punto 3"
        self.T[3]=self.T_gc+self.T_ref
        self.P[3]=self.P_gc*10**5   #iperparametro
        self.mix.update(CP.PT_INPUTS, self.P[3], self.T[3])
        self.H[3]=self.mix.hmass()
        self.S[3]=self.mix.smass()
         
          
        "punto 7"
        self.T[7]=self.T_eva+self.T_ref
        self.mix.update(CP.QT_INPUTS, 1, self.T[7])
        self.P[7]=self.mix.p()
        self.H[7]=self.mix.hmass()
        H7_tot=self.H[7]+1/2*self.v**2
        D7=self.mix.rhomass()
           
        "punto 0"
        self.T[0]=self.T_EVA_AT+self.T_ref
        self.mix.update(CP.QT_INPUTS, 1, self.T[0])
        self.P[0]=self.mix.p()
        self.H[0]=self.mix.hmass()
        
        "punto 1"
        self.P[1]=self.P[0]
        self.T[1]=self.T[0] + self.eps*(self.T[3]-self.T[0])
        self.mix.update(CP.PT_INPUTS, self.P[1], self.T[1])
        self.H[1]=self.mix.hmass()
        self.S[1]=self.mix.smass()
        
        "PUNTO 2"
        self.P[2]=self.P[3]
        self.mix.update(CP.PSmass_INPUTS, self.P[2], self.S[1])
        H2_id=self.mix.hmass()
        self.H[2]=self.H[1]+(H2_id-self.H[1])/self.eta_c
        
        "punto 4"
        self.P[4]=self.P[3]
        self.H[4]=self.H[3]-self.H[1]+self.H[0]
        self.mix.update(CP.HmassP_INPUTS, self.H[4], self.P[4])
        self.S[4]=self.mix.smass()
        
        "punto 9"
        self.P[9]=self.P[0]
        
        'primario'
        self.P[5]=self.P[7]-1/2*D7*self.v**2
        self.mix.update(CP.PSmass_INPUTS, self.P[5], self.S[4])
        H5_is=self.mix.hmass()
        
        self.H[5]=self.H[4]-self.eta1*(self.H[4]-H5_is)
        
        H4_tot=self.H[4]+1/2*self.v**2
        v_5=(2*(H4_tot-self.H[5]))**(1/2)
        
        
        
        "IPOTIZZO"
        
        def loop(omega):
            'mixing'
            v_mix=self.eta_mix*(v_5+ omega*self.v)/(1+omega)
            H8_tot=(H4_tot+omega*H7_tot)/(1+omega)
            #H8_tot=H4_tot+omega*H7_tot sbagliato
            self.H[8]=H8_tot - 1/2*v_mix**2
            self.P[8]=self.P[5]
            self.mix.update(CP.HmassP_INPUTS, self.H[8], self.P[8])
            self.S[8]=self.mix.smass()
            
            
            'diffusore'
            H9_tot=H8_tot
            self.H[9]=H9_tot - 1/2*self.v**2
            
            H9_is= self.H[8] + self.eta2*(self.H[9]-self.H[8])
# =============================================================================
#                 self.mix.update(CP.HmassSmass_INPUTS, H9_is, self.S[8])
#                 P9_new=mix.p()
# =============================================================================
            "diffusore, ragioniamo all'incontrario per evitare HmassSmass_INPUTS"
            self.mix.update(CP.HmassP_INPUTS, H9_is, self.P[9])
            S9_is=self.mix.smass()
        
            #print(P[10]*10**-5,'\n',Q) 
            return S9_is-self.S[8]
                
        omega=fsolve(loop,0.5)
        self.m[1]=1/(1+omega)
        self.m[2]=omega*self.m[1]
        
        "PUNTI RIMANENTI"
        
        "punto 6"
        self.P[6]=self.P[7]
        self.H[6]=self.H[4]
        self.mix.update(CP.HmassP_INPUTS, self.H[6], self.P[6])
        self.T[6]=self.mix.T()
       
        
        "punto 8"
        self.mix.update(CP.HmassP_INPUTS, self.H[8], self.P[8])
        self.T[8]=self.mix.T()
        
        "punto 5"
        self.mix.update(CP.HmassP_INPUTS, self.H[5], self.P[5])
        self.T[5]=self.mix.T()
        
        "punto 9"
        self.mix.update(CP.HmassP_INPUTS,self.H[9], self.P[9])
        self.T[9]=self.mix.T()
        
        
        cop=(self.H[1]-self.H[3])/(self.H[2]-self.H[1])
        
        
    
        return self.T,self.P,self.H,self.S,self.m,cop



    def punti_rimanenti(self,T,P,H,S):    
        
        
        "PUNTO 6"
        T[6]=T[10]
        P[6]=P[10]
        self.mix_l.update(CP.PQ_INPUTS, P[6], 0)
        self.mix_l.Q()
        H[6]=self.mix_l.hmass()
        
        
        "PUNTO 7"
        H[7]=H[6]
        P[7]=P[8]
        self.mix_l.update(CP.HmassP_INPUTS, H[7], P[7])
        T[7]=self.mix_l.T()
        #S[7]=mix_l.smass()
        
        
        
        return T,P,H,S
    
    def entropie(self,T,P,H,S):    
        "PUNTO 4"
        self.mix_g.update(CP.HmassP_INPUTS, H[4], P[4])
        T[4]=self.mix.T()
        
        "punto 5"
        self.mix_g.update(CP.PSmass_INPUTS, P[5], S[4])
        T[5]=self.mix_g.T()

        self.mix_g.update(CP.HmassP_INPUTS, H[2], P[2])
        T[2]=self.mix_g.T()
        S[2]=self.mix_g.smass()
                
        return T,P,H,S

if __name__ == "__main__":
    
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
    
    x= 0.96
    #mix.set_mole_fractions([x, 1-x])
    mix.set_mass_fractions([x, 1-x])
    
    
    eps=0.8
    T_gc=40
    P_gc=90
    T_eva=-5
    T_EVA_AT=T_eva+3#7
    eta_c=0.8
    eta1=0.9
    eta2=0.8
    eta_mix=0.85
    v=10 #m/s

    T_ref=273.15
    

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
    plt.plot((0,(H[0]-H[9]+(H[7]-H[6])*m[2])/1000),(T[6]+5-T_ref,T[0]+5-T_ref),'r')
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("T [°C]")
    plt.title('Evaporatore')
    plt.grid()

