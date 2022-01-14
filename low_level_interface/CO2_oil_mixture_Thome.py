import numpy as np 
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import compressore as c
import grafici_termodinamici as gt



a=[-2394.5, 182.52, -724.21, 3868.0, -5268.9]
b=[8.0736, -0.72212, 2.3914, -13.779, 17.066]

#P_sat=0.55#la vuole in MPa

def T_bub(w,P_sat):
    P_sat=P_sat/1000000
    
    A=a[0]+a[1]*w+a[2]*w**3+a[3]*w**5+a[4]*w**7
    B=b[0]+b[1]*w+b[2]*w**3+b[3]*w**5+b[4]*w**7

    return A/(np.log(P_sat)-B)


    
T=np.zeros(4)
P=np.zeros(4)
H=np.zeros(4)
H_f=np.zeros(4)
H_o=np.zeros(4)
Q=np.zeros(4)
S=np.zeros(4)
w=np.ones(4)


T_ref=273.15


lib="REFPROP"
fld   = CP.AbstractState(lib, "CO2")  
oil   = CP.AbstractState("INCOMP", "SAB")


"punto 2"  
k=0.03 #fraz di olio nel totale   massa olio/massa totale
#w =frazione di olio nel liquido  massa olio/massa liquida 
T[2]=40+T_ref
P[2]=95*10**5
fld.update(CP.PT_INPUTS, P[2], T[2])
H_f[2]=fld.hmass()

oil.update(CP.PT_INPUTS, P[2], T[2])
H_o[2]=oil.hmass()

H[2]=k*H_o[2] + (1-k)*H_f[2]


"punto 3"
H[3]=H[2]
T[3]=-5+T_ref

fld.update(CP.HmassT_INPUTS, H[3], T[3])
P[3]=fld.p()    # P di tentativo

P1=P[3]+100
P2= P[3]-100
fld.update(CP.PQ_INPUTS, P1 , 0)
T1= fld.T()

fld.update(CP.PQ_INPUTS, P2 , 0)
T2= fld.T()

"trovo a0 e b0 considerando w=0"
P1=P1*10**-6
P2=P2*10**-6

b[0]=(T2*np.log(P2)-T1*np.log(P1))/(T2-T1)

a[0]=T1*(np.log(P1)-b[0])


def P_3(p3):
    P[3]=p3
    
    oil.update(CP.PT_INPUTS, P[3], T[3])
    H_o[3]=oil.hmass()

    H_f[3]=( H[3]- k*H_o[3] )/(1-k)
    
    fld.update(CP.HmassT_INPUTS, H_f[3], T[3])
    Q[3]=fld.Q()


    z_l=(1-k-Q[3]+k*Q[3])/(1-Q[3]+k*Q[3])

    w[3]=1-z_l

    #print(p3)

    return T_bub(w[3],P[3])-T[3]

P[3]=fsolve(P_3, P[3])


"punto 0"
P[0]=P[3]
w[0]= 0.56  #è ottimizzabile siamo noi a deciderlo
Q_tot=1-k/w[0]
z_l=1-w[0]
Q[0]=(1-z_l-k)/(1-z_l-k+z_l*k)



T[0]=T_bub(w[0],P[0])


oil.update(CP.PT_INPUTS, P[0], T[0])
H_o[0]=oil.hmass()

fld.update(CP.QT_INPUTS, Q[0], T[0])
H_f[0]=fld.hmass()

H[0] = k*H_o[0] + (1-k)*H_f[0]


"punto 1"
P[1]=P[2]
T[1]=c.Temperatura_mandata(T[3]-T_ref, T[0]-T_ref, P[1]*10**-5)+T_ref

fld.update(CP.PT_INPUTS, P[1], T[1])
H_f[1]=fld.hmass()

oil.update(CP.PT_INPUTS, P[1], T[1])
H_o[1]=oil.hmass()

H[1]=k*H_o[1] + (1-k)*H_f[1]


COP=(H[0]-H[2])/(H[1]-H[0])
print('cop con olio='+str(COP))



"CICLO SOLO CO2"
t=np.zeros(4)
p=np.zeros(4)
h=np.zeros(4)

"punto 2"  
t[2]=T[2]
p[2]=P[2]
fld.update(CP.PT_INPUTS, p[2], t[2])
h[2]=fld.hmass()

"punto 3"
h[3]=h[2]
t[3]=T[3]
fld.update(CP.HmassT_INPUTS, h[3], t[3])
p[3]=fld.p()    # P di tentativo

"punto 0"
p[0]=p[3]
fld.update(CP.PQ_INPUTS, p[0], 1)
t[0]=fld.T()
h[0]=fld.hmass()

"punto 1"
p[1]=p[2]
t[1]=c.Temperatura_mandata(t[3]-T_ref, t[0]-T_ref, p[1]*10**-5)+T_ref

fld.update(CP.PT_INPUTS, p[1], t[1])
h[1]=fld.hmass()

cop=(h[0]-h[2])/(h[1]-h[0])
print('cop senza olio='+str(cop))


plt.figure(dpi=200)
gt.grafico_PH_super_semplice(P/100000,H/1000, 'r',0)
gt.grafico_PH_super_semplice(p/100000,h/1000, 'b',1)