import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c
import matplotlib.pyplot as plt


class Funz:  
    def __init__(self,eps,P_gc,T_gc,T_eva,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g):        
        

        self.mix = mix
        self.mix_l = mix_l
        self.mix_g = mix_g

        self.eps=eps
        self.T_gc=T_gc
        self.P_gc=P_gc
        self.T_eva=T_eva
        self.eta_c=eta_c
        self.eta1=eta1
        self.eta2=eta2
        self.eta_mix=eta_mix
        self.v=v
       
        self.T=np.zeros(11)
        self.P=np.zeros(11)
        self.H=np.zeros(11)
        self.S=np.zeros(11)
        self.m=np.ones(3)
        
        
        self.T_ref=273.15
               


    def imp(self):
        
        "PUNTI NOTI A PRIORI" #non ci sono
        "punto 3"
        self.T[3]=self.T_gc+self.T_ref
        self.P[3]=self.P_gc*10**5   #iperparametro
         
          
        "punto 8"
        self.T[8]=self.T_eva+self.T_ref
        
        "punto 4"
        self.P[4]=self.P[3]
        
        "punto 2"
        self.P[2]=self.P[3]
            
        def loop(x):
            #print(x)
            f=np.zeros(2)
                
            
            "VARIABILI"
            self.P[10]=x[0]*10**5 #3642449.4142085267#
            Q=x[1] #titolo molare in 10 #0.8252851389039894#
            
            "punto 10"
            self.mix.update(CP.PQ_INPUTS, self.P[10], Q)
            self.T[10]=self.mix.T()
            H10_in=self.mix.hmass()
            self.S[10]=self.mix.smass()
            
            "SEPARATORE"
            m_l=self.mix.mole_fractions_liquid()
            m_g=self.mix.mole_fractions_vapor()
           
            
            self.mix_l.set_mole_fractions(m_l)
            self.mix_g.set_mole_fractions(m_g)
            
            "punto 3"
            self.mix_g.update(CP.PT_INPUTS, self.P[3], self.T[3])
            self.H[3]=self.mix_g.hmass()
            
            "PUNTO 8"
            self.mix_l.update(CP.QT_INPUTS, 1, self.T[8])
            self.P[8]=self.mix_l.p()
            self.H[8]=self.mix_l.hmass()
            D8=self.mix_l.rhomass()
            H8_tot=self.H[8]+1/2*self.v**2
            
            "PUNTO 0"
            self.T[0]=self.T[10]
            self.P[0]=self.P[10]
            self.mix_g.update(CP.PQ_INPUTS, self.P[0], 1)
            #mix_g.update(CP.PT_INPUTS, P[0], T[0])
            #self.mix_g.Q()
            self.H[0]=self.mix_g.hmass()
            
            "punto 1"
            self.P[1]=self.P[0]
            self.T[1]=self.T[0] + self.eps*(self.T[3]-self.T[0])
            self.mix_g.update(CP.PT_INPUTS, self.P[1], self.T[1])
            self.H[1]=self.mix_g.hmass()
            
            "PUNTO 4"
            self.H[4]=self.H[3]-self.H[1]+self.H[0]
            self.mix_g.update(CP.HmassP_INPUTS, self.H[4], self.P[4])
            self.S[4]=self.mix_g.smass()
            
            'primario'
            self.P[5]=self.P[8]-1/2*D8*self.v**2
            self.mix_g.update(CP.PSmass_INPUTS, self.P[5], self.S[4])
            H5_is=self.mix_g.hmass()
            
            self.H[5]=self.H[4]-self.eta1*(self.H[4]-H5_is)
            
            H4_tot=self.H[4]+1/2*self.v**2
            v_5=(2*(H4_tot-self.H[5]))**(1/2)
            
            "CIRCUITO INTERNO"
            q=Q*self.mix_g.molar_mass()/self.mix.molar_mass() #massico
            omega=1/q-1
            self.m[1]=1-q
            self.m[2]=q
                
            'mixing'
            v_mix=self.eta_mix*(v_5+ omega*self.v)/(1+omega)
            H9_tot=(H4_tot+omega*H8_tot)/(1+omega)
            self.H[9]=H9_tot - 1/2*v_mix**2
            self.P[9]=self.P[5] 
            self.mix.update(CP.HmassP_INPUTS, self.H[9], self.P[9])
            self.S[9]=self.mix.smass()
            
            
            'diffusore'
            H10_tot=H9_tot
            self.H[10]=H10_tot - 1/2*self.v**2
            
            H10_is= self.H[9] + self.eta2*(self.H[10]-self.H[9])
        # =============================================================================
        #     #print('aaaaa')
        #     mix.update(CP.HmassSmass_INPUTS, H10_is, S[9])
        #     #print('bbbbb')
        #     
        #     
        #     
        #     
        # 
        #     #mix.update(CP.HmassSmass_INPUTS, H[10], S[9])
        #     P[10]=mix.p()
        #     Q=mix.Q()
        # =============================================================================
                
            "diffusore, ragioniamo all'incontrario per evitare HmassSmass_INPUTS"
            self.mix.update(CP.HmassP_INPUTS, H10_is, self.P[10])
            S10_is=self.mix.smass()
        
            #print(P[10]*10**-5,'\n',Q) 
            f[0]=S10_is-self.S[9]
            f[1]=self.H[10]-H10_in
            #f[0]=P[10]-x[0]*10**5
            #f[1]=Q-x[1]
            
            return f
        
        fsolve(loop,[30,0.7])

  
        
        "punto 1"
        self.mix_g.update(CP.PT_INPUTS, self.P[1], self.T[1])
        self.S[1]=self.mix_g.smass()
        
        "PUNTO 2"
        self.mix_g.update(CP.PSmass_INPUTS, self.P[2], self.S[1])
        H2_id=self.mix_g.hmass()
        self.H[2]=self.H[1]+(H2_id-self.H[1])/self.eta_c
        
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
        self.mix_g.update(CP.PSmass_INPUTS, P[5], S[4])#mi sa sbagliato casomai H5
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
    
    x= 1#0.979#0.95#0.8904072#0.890685902909111#0.8906859#0.94470172#0.99
    #mix.set_mole_fractions([x, 1-x])
    mix.set_mass_fractions([x, 1-x])
    
    
    eps=0.8
    T_gc=40
    P_gc=90
    T_eva=-5-2.980773498249846#-5
    eta_c=0.8
    eta1=0.9
    eta2=0.8
    eta_mix=0.85
    v=10 #m/s

    T_ref=273.15
    

    funz=Funz(eps,P_gc,T_gc,T_eva,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    
    print(cop)
    
    
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    T,P,H,S=funz.entropie(T,P,H,S)

    cop2=m[1]*(H[8]-H[7])/(H[2]-H[1])/m[2]
    print(cop2)
    print(T[8]-T[7])
# =============================================================================
#     pc=60*10**5
#     gt.grafico_PH_eiettore_sep(P/100000,H/1000, mix, mix_l, mix_g,pc, "CO2&R1234YF",x)
# 
#    
#     "evaporatore"
#     n=100
#     H_e=np.linspace(H[7],H[8],n)
#     T_e=np.zeros(n)
#     for i in range(n):
#         mix_l.update(CP.HmassP_INPUTS, H_e[i], P[7])
#         T_e[i]=mix_l.T()
#     
#     plt.figure(dpi=200)
#     plt.plot(H_e/1000,T_e-T_ref)
#     plt.plot((H[7]/1000,H[8]/1000),(T[7]+5-T_ref,T[8]+5-T_ref),'r')
#     plt.xlabel("H [kJ/kg]")
#     plt.ylabel("T [°C]")
#     plt.title('Evaporatore')
#     plt.grid()
# 
# =============================================================================
