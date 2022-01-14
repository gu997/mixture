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

x=0.9
#mix.set_mole_fractions([x, 1-x])
mix.set_mass_fractions([x, 1-x])


eps=0.8
T_gc=40
P_gc=95
T_eva=-5
T_sep=10
eta_c=0.8


T=np.zeros(11)
P=np.zeros(11)
H=np.zeros(11)
S=np.zeros(11)
m=np.ones(3)
m_ll=np.zeros(2)
m_gg=np.zeros(2)


T_ref=273.15


"PUNTI NOTI A PRIORI"

"punto 3"
T[3]=T_gc+T_ref
P[3]=P_gc*10**5   #iperparametro
mix.update(CP.PT_INPUTS, P[3], T[3])
H[3]=mix.hmass()

"punto 8"
T[8]=T_eva+T_ref
        
"punto 5"
T[5]=T_sep+T_ref

"punto 4"
P[4]=P[3]

"punto 2"
P[2]=P[3]

"IPOTIZZO IL titolo in  5"
def loop(Q):
    
    "punto 5"
    #Q=0.3 #molare
    mix.update(CP.QT_INPUTS, Q, T[5])
    P[5]=mix.p()
    H[5]=mix.hmass()
    
    "SEPARATORE"
    m_l=mix.mole_fractions_liquid()
    m_g=mix.mole_fractions_vapor()
    m_ll[0]=m_l[0]
    m_gg[0]=m_g[0]
    
    mix_l.set_mole_fractions(m_l)
    mix_g.set_mole_fractions(m_g)
    
    "PUNTO 8"
    mix_l.update(CP.QT_INPUTS, 1, T[8])
    P[8]=mix_l.p()
    H[8]=mix_l.hmass()
    
    "PUNTO 9"
    T[9]=T[5]
    P[9]=P[5]
    mix_g.update(CP.PT_INPUTS, P[9], T[9])
    mix_g.Q()
    H[9]=mix_g.hmass()
    
    "PUNTO 10"
    H[10]=H[9]
    P[10]=P[8]
    #mix_l.update(CP.HmassP_INPUTS, H[10], P[10])
    #T[10]=mix_g.T()
    
    "PUNTO 0"
    q=Q*mix_g.molar_mass()/mix.molar_mass() #massico
    m[1]=1-q
    m[2]=q
    H[0]=m[1]*H[8]+m[2]*H[10]
    P[0]=P[8]
    mix.update(CP.HmassP_INPUTS, H[0], P[0])
    T[0]=mix.T()
    #S[0]=mix.smass()
    
    "punto 1"
    P[1]=P[8]
    T[1]=T[0] + eps*(T[3]-T[0])
    mix.update(CP.PT_INPUTS, P[1], T[1])
    H[1]=mix.hmass()
    
    "PUNTO 4"
    H[4]=H[3]-H[1]+H[0]
    
    "punto 5"
    H[5]=H[4]
    mix.update(CP.HmassT_INPUTS, H[5], T[5])

   
    return mix.Q()-Q

Q=fsolve(loop,0.3)

"PUNTO 6"
P[6]=P[5]
T[6]=T[5]
mix_l.update(CP.PT_INPUTS, P[6], T[6])
mix_l.Q()
H[6]=mix_l.hmass()
#S[6]=mix_l.smass()


"PUNTO 7"
H[7]=H[6]
P[7]=P[8]
mix_l.update(CP.HmassP_INPUTS, H[7], P[7])
T[7]=mix_l.T()
#S[7]=mix_l.smass()

"PUNTO 4"
mix.update(CP.HmassP_INPUTS, H[4], P[4])
T[4]=mix.T()


"PUNTO 10"
mix_g.update(CP.HmassP_INPUTS, H[10], P[10])
T[10]=mix_g.T()

# =============================================================================
# "punto 2"
# T[2]=c.Temperatura_mandata(T[0]-T_ref, T[1]-T_ref, P[2]*10**-5)+T_ref
# mix.update(CP.PT_INPUTS, P[2], T[2])
# H[2]=mix.hmass()
# =============================================================================


"punto 1"
mix.update(CP.PT_INPUTS, P[1], T[1])
S[1]=mix.smass()

"PUNTO 2"
mix.update(CP.PSmass_INPUTS, P[2], S[1])
H2_id=mix.hmass()
H[2]=H[1]+(H2_id-H[1])/eta_c
mix.update(CP.HmassP_INPUTS, H[2], P[2])
T[2]=mix.T()
S[2]=mix.smass()



cop=m[1]*(H[8]-H[7])/(H[2]-H[1])
cop2=(H[1]-H[3])/(H[2]-H[1])
print('cop =',cop)
print('cop2=',cop2)

pc=65*10**5

gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g,pc, "CO2&R1234YF",x)
