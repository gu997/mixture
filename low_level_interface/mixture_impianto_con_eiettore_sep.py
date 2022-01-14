import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c



lib="REFPROP"

mix   = CP.AbstractState(lib, "CO2&R1234YF")
mix_l = CP.AbstractState(lib, "CO2&R1234YF")
mix_g = CP.AbstractState(lib, "CO2&R1234YF")

# =============================================================================
# mix   = CP.AbstractState(lib, "CO2&CO2")
# mix_l = CP.AbstractState(lib, "CO2&CO2")
# mix_g = CP.AbstractState(lib, "CO2&CO2")
# =============================================================================

x=1#0.99#0.99#0.9
mix.set_mole_fractions([x, 1-x])
#mix.set_mass_fractions([x, 1-x])


eps=0.8
T_gc=40
P_gc=90
T_eva=-8.645466228845635#-5
eta_c=0.8
eta1=0.9
eta2=0.8
eta_mix=0.85
v=10 #m/s

T=np.zeros(11)
P=np.zeros(11)
H=np.zeros(11)
S=np.zeros(11)
m=np.ones(3)
m_ll=np.zeros(2)
m_gg=np.zeros(2)


T_ref=273.15


"PUNTI NOTI A PRIORI" #non ci sono
"punto 3"
T[3]=T_gc+T_ref
P[3]=P_gc*10**5   #iperparametro
 
  
"punto 8"
T[8]=T_eva+T_ref

"punto 4"
P[4]=P[3]

"punto 2"
P[2]=P[3]
    
def loop(x):
    print(x)
    f=np.zeros(2)
        
    
    "VARIABILI"
    P[10]=x[0]*10**5 #3642449.4142085267#
    Q=x[1] #titolo molare in 10 #0.8252851389039894#
    
    "punto 10"
    mix.update(CP.PQ_INPUTS, P[10], Q)
    T[10]=mix.T()
    H10_in=mix.hmass()
    S[10]=mix.smass()
    
    "SEPARATORE"
    m_l=mix.mole_fractions_liquid()
    m_g=mix.mole_fractions_vapor()
    m_ll[0]=m_l[0]
    m_gg[0]=m_g[0]
    
    mix_l.set_mole_fractions(m_l)
    mix_g.set_mole_fractions(m_g)
    
    "punto 3"
    mix_g.update(CP.PT_INPUTS, P[3], T[3])
    H[3]=mix_g.hmass()
    
    "PUNTO 8"
    mix_l.update(CP.QT_INPUTS, 1, T[8])
    P[8]=mix_l.p()
    H[8]=mix_l.hmass()
    D8=mix_l.rhomass()
    H8_tot=H[8]+1/2*v**2
    
    "PUNTO 0"
    T[0]=T[10]
    P[0]=P[10]
    mix_g.update(CP.PQ_INPUTS, P[0], 1)
    #mix_g.update(CP.PT_INPUTS, P[0], T[0])
    mix_g.Q()
    H[0]=mix_g.hmass()
    
    "punto 1"
    P[1]=P[0]
    T[1]=T[0] + eps*(T[3]-T[0])
    mix_g.update(CP.PT_INPUTS, P[1], T[1])
    H[1]=mix_g.hmass()
    
    "PUNTO 4"
    H[4]=H[3]-H[1]+H[0]
    mix_g.update(CP.HmassP_INPUTS, H[4], P[4])
    S[4]=mix_g.smass()
    
    'primario'
    P[5]=P[8]-1/2*D8*v**2
    mix_g.update(CP.PSmass_INPUTS, P[5], S[4])
    H5_is=mix_g.hmass()
    
    H[5]=H[4]-eta1*(H[4]-H5_is)
    
    H4_tot=H[4]+1/2*v**2
    v_5=(2*(H4_tot-H[5]))**(1/2)
    
    "CIRCUITO INTERNO"
    q=Q*mix_g.molar_mass()/mix.molar_mass() #massico
    omega=1/q-1
    m[1]=1-q
    m[2]=q
        
    'mixing'
    v_mix=eta_mix*(v_5+ omega*v)/(1+omega)
    H9_tot=(H4_tot+omega*H8_tot)/(1+omega)
    H[9]=H9_tot - 1/2*v_mix**2
    P[9]=P[5] 
    mix.update(CP.HmassP_INPUTS, H[9], P[9])
    S[9]=mix.smass()
    
    
    'diffusore'
    H10_tot=H9_tot
    H[10]=H10_tot - 1/2*v**2
    
    H10_is= H[9] + eta2*(H[10]-H[9])
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
    mix.update(CP.HmassP_INPUTS, H10_is, P[10])
    S10_is=mix.smass()

    #print(P[10]*10**-5,'\n',Q) 
    f[0]=S10_is-S[9]
    f[1]=H[10]-H10_in
    #f[0]=P[10]-x[0]*10**5
    #f[1]=Q-x[1]
    
    return f

fsolve(loop,[30,0.7])


"PUNTO 6"
T[6]=T[10]
P[6]=P[10]
mix_l.update(CP.PQ_INPUTS, P[6], 0)
mix_l.Q()
H[6]=mix_l.hmass()


"PUNTO 7"
H[7]=H[6]
P[7]=P[8]
mix_l.update(CP.HmassP_INPUTS, H[7], P[7])
T[7]=mix_l.T()
#S[7]=mix_l.smass()

"PUNTO 4"
mix_g.update(CP.HmassP_INPUTS, H[4], P[4])
T[4]=mix.T()

"punto 5"
mix_g.update(CP.PSmass_INPUTS, P[5], S[4])
T[5]=mix_g.T()

"punto 9"
mix.update(CP.HmassP_INPUTS, H[9], P[9])
T[9]=mix_g.T()

"punto 1"
mix_g.update(CP.PT_INPUTS, P[1], T[1])
S[1]=mix_g.smass()

"PUNTO 2"
mix_g.update(CP.PSmass_INPUTS, P[2], S[1])
H2_id=mix_g.hmass()
H[2]=H[1]+(H2_id-H[1])/eta_c
mix_g.update(CP.HmassP_INPUTS, H[2], P[2])
T[2]=mix_g.T()
S[2]=mix_g.smass()



cop=m[1]*(H[8]-H[7])/(H[2]-H[1])/m[2]
cop2=(H[1]-H[3])/(H[2]-H[1])
print('cop =',cop)
print('cop2=',cop2)

pc=65*10**5

#gt.grafico_PH_eiettore_sep(P/100000,H/1000, mix, mix_l, mix_g,pc, "CO2&R1234YF",x)
