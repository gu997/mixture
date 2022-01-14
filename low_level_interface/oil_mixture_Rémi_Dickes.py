import CoolProp.CoolProp as CP
import numpy as np 
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


"R245fa-oil mixture"

a1 = -5.47134e-02
a2 =  5.62233e-02
a3 =  2.36502e-02
a4 = -4.47031e-02
a5 =  2.91990e-02
a6 = -7.20697e-03
a7 =  6.32951e-04

#lib="HEOS"
lib="REFPROP"
fld   = CP.AbstractState(lib, "R245fa")
oil   = CP.AbstractState("INCOMP", "SAB")

eta_c=0.2

# =============================================================================
# m_f_liq=1
# m_oil=0.03
# 
# z_l= m_f_liq/(m_f_liq + m_oil)  #frazione di fliudo di lavoro nella fase liquida
# =============================================================================

def h_oil(T, P):
# =============================================================================
#     a2_oil = 699.4
#     a3_oil = 3.976
#      
#      
#      
#     #a2_oil = 1158.8
#     #a3_oil = 2.3639
#     h=a2_oil*T + a3_oil/2*T**2
# =============================================================================
    
    oil.update(CP.PT_INPUTS, P , T)
    
    h=oil.hmass()
    
    return h

def S_oil(T, P):
    oil.update(CP.PT_INPUTS, P , T)
    
    s=oil.smass()
    
    return s

def t( p, z_l):
    
    A = a1 + a2/(z_l**(1/2))
    B=a3 + a4/(z_l**(1/2)) + a5/z_l + a6/(z_l**(3/2)) + a7/(z_l**2)

    teta = (1-z_l)*(A+B*p*10**-5)
    
    fld.update(CP.PQ_INPUTS, p, 0)
    Tsat=fld.T()
    tt=Tsat*(1+teta)
    #print(tt)
    return tt
    
 
def t_0(p , args):
    return t( p, z_l)-args



"teta = (T − Tsat(P) )/Tsat(P)"
"allora  T = Tsat(P)* (1+teta)"

"h_mix=mf_l*Hf_l + mf_v*Hf_v + m_oil*H_oil"
"h_mix= (1-Q)*(1-k)*Hf_l + Q*(1-k)*Hf_v + k*H_oil"
"Q=mf_v/mf"


T=np.zeros(4)
P=np.zeros(4)
H=np.zeros(4)
Q=np.zeros(4)
S=np.zeros(4)




T_ref=273.15

k=0.03 #frazione di olio nel totale

"punto 2"
T[2]=40+T_ref
#print('T2  ',T[2])
z_l = 1-k #poichè tutto il fluido è liquido
fld.update(CP.QT_INPUTS, 0, T[2])
P[2]=fld.p()    # P di tentativo
H[2] = k*h_oil(T[2], P[2]) + (1-k)*fld.hmass()

P[2] = fsolve(t_0, [P[2]],args=(T[2]))
#t(P[2],z_l)

"punto 3"
T[3]=0+T_ref
H[3]=H[2]



fld.update(CP.QT_INPUTS, 0, T[3])
P[3]=fld.p()    # P di tentativo

H_f= H[3]- k*h_oil(T[3], P[3])

H_f_l = fld.hmass()


fld.update(CP.QT_INPUTS, 1, T[3])
H_f_v = fld.hmass()

Q[3] = ( H_f - (1-k)*H_f_l )/(1-k)/(H_f_v-H_f_l)

z_l=(1-k-Q[3]+k*Q[3])/(1-Q[3]+k*Q[3])

P[3] = fsolve(t_0, [P[3]],args=(T[3]))


"punto 0"
P[0]=P[3]
Q[0]=0.9 #○altrimenti z_l vale 0

z_l=(1-k-Q[0]+k*Q[0])/(1-Q[0]+k*Q[0])
T[0]=t(P[0],z_l)


fld.update(CP.QT_INPUTS, 0, T[0])
H_f_l = fld.hmass()
S_f_l = fld.smass()

fld.update(CP.QT_INPUTS, 1, T[0])
H_f_v = fld.hmass()
S_f_v = fld.smass()

H[0]= (1-Q[0])*(1-k)*H_f_l + Q[0]*(1-k)*H_f_v + k*h_oil(T[0], P[0])
H0  = (1-k)*H_f_v + k*h_oil(T[0], P[0])  #faccio finta che sia tutto vapore

S[0]= (1-Q[0])*(1-k)*S_f_l + Q[0]*(1-k)*S_f_v + k*S_oil(T[0], P[0])
print(fld.Q())

"punto 1"

P[1]=P[2]
S[1]=S[0]
#print('S0=',S[0])
#metto T di tentativo
fld.update(CP.PSmass_INPUTS, P[1], S[1])
T[1]=fld.T()

def isos(Tt):
    #print(Tt)
    S_f=S[1]- k*S_oil(Tt, P[1])
    #print('S_f=',S_f)
    fld.update(CP.PSmass_INPUTS, P[1], S_f)
    #print('fld.T()=',fld.T())
    return fld.T()-Tt

T[1]=fsolve(isos,T[1])
H1_id=fld.hmass()

H[1]=H[0]+(H1_id-H[0])/eta_c
fld.update(CP.HmassP_INPUTS, H[1], P[1])
T[1]=fld.T()
print(fld.Q())

plt.figure(dpi=200)
plt.plot(H/1000,P/100000)

plt.figure(dpi=200)
plt.plot(H/1000,T-T_ref)