import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c

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
        
        "punto 7"
        self.T[7]=self.T_eva+self.T_ref
                
        "punto 5"
        self.T[5]=self.T_sep+self.T_ref
        
        "punto 4"
        #self.P[4]=self.P[3]
        
        "punto 2"
        self.P[2]=self.P[3]
        
        "IPOTIZZO IL titolo in  5"
        def loop(Q):
            
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
            
            "PUNTO 6"
            self.P[6]=self.P[5]
            self.T[6]=self.T[5]
            self.mix_l.update(CP.PT_INPUTS, self.P[6], self.T[6])
            #mix_l.Q()
            self.H[6]=self.mix_l.hmass()
            #S[6]=mix_l.smass()
        
            "PUNTO 7"
            self.H[7]=self.H[6]
            self.mix_l.update(CP.HmassT_INPUTS, self.H[7], self.T[7])
            self.P[7]=self.mix_l.p()
            
            "PUNTO 8"
            self.P[8]=self.P[7]
            self.mix_l.update(CP.PQ_INPUTS,  self.P[8], 1)
            #self.P[8]=self.mix_l.p()
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
# =============================================================================
#         "PUNTO 6"
#         P[6]=P[5]
#         T[6]=T[5]
#         self.mix_l.update(CP.PT_INPUTS, P[6], T[6])
#         #mix_l.Q()
#         H[6]=self.mix_l.hmass()
#         #S[6]=mix_l.smass()
#         
#         
#         "PUNTO 7"
#         H[7]=H[6]
#         P[7]=P[8]
#         self.mix_l.update(CP.HmassP_INPUTS, H[7], P[7])
#         T[7]=self.mix_l.T()
#         #S[7]=mix_l.smass()
# =============================================================================

        "PUNTO 8"
        self.mix_l.update(CP.PQ_INPUTS,  self.P[8], 1)
        self.T[8]=self.mix_l.T()
        
        "PUNTO 4"
        P[4]=P[3]
        self.mix.update(CP.HmassP_INPUTS, H[4], P[4])
        T[4]=self.mix.T()
        
        
        "PUNTO 10"
        P[10]=P[8]
        self.mix_g.update(CP.HmassP_INPUTS, H[10], P[10])
        T[10]=self.mix_g.T()
        
        return T,P,H,S


if __name__ == "__main__":
    
    lib="REFPROP"

    fluids="CO2&R1234YF"
    #fluids="CO2&PROPANE"
    #fluids="CO2&ISOBUTANE"
    #fluids="CO2&HEXANE"
    


    mix   = CP.AbstractState(lib, fluids)
    mix_l = CP.AbstractState(lib, fluids)
    mix_g = CP.AbstractState(lib, fluids)
    
    x=0.9
    #mix.set_mole_fractions([x, 1-x])
    mix.set_mass_fractions([x, 1-x])
    
    
    eps=0.8
    T_gc=40
    P_gc=95
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

    pc=60*10**5
    gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g,pc, fluids,x)
    #gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g)

