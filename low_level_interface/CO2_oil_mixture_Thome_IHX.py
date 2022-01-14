"non considero separatore olio"

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

def f_T_bub(P_sat,w):
    
    return T_bub(w,P_sat)-T[0]

# =============================================================================
# def h_oil(T): #preso da Rémi Dickes
#     a2_oil = 699.4
#     a3_oil = 3.976
#      
#      
#      
#     #a2_oil = 1158.8
#     #a3_oil = 2.3639
#     h=a2_oil*T + a3_oil/2*T**2
#     
#     return h
# =============================================================================

def h_oil(T): #preso da Rémi Dickes
    a2_oil = 699.4
    a3_oil = 3.976
     
     
     
    #a2_oil = 1158.8
    #a3_oil = 2.3639
    h_r=a2_oil*273.15 + a3_oil/2*273.15**2
    h=a2_oil*T + a3_oil/2*T**2  -h_r
    
    return h/1000 + 200000



    
T=np.zeros(6)
P=np.zeros(6)
H=np.zeros(6)
H_f=np.zeros(6)
H_o=np.zeros(6)
Q=np.zeros(6)
S=np.zeros(6)
w=np.ones(6)

eps=0.8

T_ref=273.15


lib="REFPROP"
fld   = CP.AbstractState(lib, "CO2")  
#oil   = CP.AbstractState("INCOMP", "SAB")


"punto 3"  
k=0.03 #fraz di olio nel totale   massa olio/massa totale
#w =frazione di olio nel liquido  massa olio/massa liquida 
T[3]=40+T_ref
P[3]=95*10**5
fld.update(CP.PT_INPUTS, P[3], T[3])
H_f[3]=fld.hmass()

#oil.update(CP.PT_INPUTS, P[3], T[3])
H_o[3]=h_oil(T[3])#oil.hmass()

H[3]=k*H_o[3] + (1-k)*H_f[3]

"punto 0"
T[0]=-5+T_ref
w[0]= 0.6  #è ottimizzabile siamo noi a deciderlo
Q_tot=1-k/w[0]
z_l=1-w[0]
Q[0]=(1-z_l-k)/(1-z_l-k+z_l*k)


#oil.update(CP.PT_INPUTS, P[0], T[0])
H_o[0]=h_oil(T[0])#oil.hmass()

fld.update(CP.QT_INPUTS, Q[0], T[0])
H_f[0]=fld.hmass()
P[0]=fld.p()    # P di tentativo


H[0] = k*H_o[0] + (1-k)*H_f[0]


P1=P[0]+100
P2= P[0]-100
fld.update(CP.PQ_INPUTS, P1 , 0)
T1= fld.T()

fld.update(CP.PQ_INPUTS, P2 , 0)
T2= fld.T()

"trovo a0 e b0 considerando w=0"
P1=P1*10**-6
P2=P2*10**-6

b[0]=(T2*np.log(P2)-T1*np.log(P1))/(T2-T1)

a[0]=T1*(np.log(P1)-b[0])

P[0]=fsolve(f_T_bub,P[0],args=(w[0]))


"punto 1"
T[1]=T[0] + eps*(T[3]-T[0])
P[1]=P[0]
fld.update(CP.PT_INPUTS, P[1], T[1])
H_f[1]=fld.hmass()

#oil.update(CP.PT_INPUTS, P[1], T[1])
H_o[1]=h_oil(T[1])#oil.hmass()

H[1]=k*H_o[1] + (1-k)*H_f[1]


"punto 4"
H[4]=H[3]-H[1]+H[0]
P[4]=P[3]
fld.update(CP.HmassP_INPUTS, H[4], P[4])
T[4]=fld.T()    # T di tentativo

def T_4(t4):
    #print(t4)

    T[4]=t4
    
    #oil.update(CP.PT_INPUTS, P[4], T[4])
    H_o[4]=h_oil(T[4])#oil.hmass()

    H_f[4]=( H[4]- k*H_o[4] )/(1-k)
    
    fld.update(CP.HmassP_INPUTS, H_f[4], P[4])

    return fld.T()-t4

T[4]=fsolve(T_4, T[4])




"punto 5"
H[5]=H[4]
P[5]=P[0]

fld.update(CP.HmassP_INPUTS, H[5], P[5])
T[5]=fld.T()    # T di tentativo

def T_5(t5):
    T[5]=t5
    
    #oil.update(CP.PT_INPUTS, P[3], T[3])
    H_o[5]=h_oil(T[5])#oil.hmass()

    H_f[5]=( H[5]- k*H_o[5] )/(1-k)
    
    fld.update(CP.HmassT_INPUTS, H_f[5], T[5])
    Q[5]=fld.Q()

    z_l=(1-k-Q[5]+k*Q[5])/(1-Q[5]+k*Q[5])

    w[5]=1-z_l

    return T_bub(w[5],P[5])-T[5]

T[5]=fsolve(T_5, T[5])


"punto 2"
P[2]=P[3]
T[2]=c.Temperatura_mandata(T[5]-T_ref, T[1]-T_ref, P[2]*10**-5)+T_ref

fld.update(CP.PT_INPUTS, P[2], T[2])
H_f[2]=fld.hmass()

#oil.update(CP.PT_INPUTS, P[2], T[2])
H_o[2]=h_oil(T[2])#oil.hmass()

H[2]=k*H_o[2] + (1-k)*H_f[2]




COP=(H[1]-H[3])/(H[2]-H[1])
print('cop con olio='+str(COP))

COP=(H[0]-H[5])/(H[2]-H[1])
print('cop con olio='+str(COP))







"CICLO SOLO CO2"







t=np.zeros(6)
p=np.zeros(6)
h=np.zeros(6)


T_ref=273.15

"PUNTI NOTI A PRIORI"

"punto 3"
t[3]=T[3]
p[3]=P[3]
fld.update(CP.PT_INPUTS, p[3], t[3])
h[3]=fld.hmass()

"punto 0"
t[0]=T[0]
fld.update(CP.QT_INPUTS, 1, t[0])
h[0]=fld.hmass()
p[0]=fld.p()    
 
"punto 1"
t[1]=t[0] + eps*(t[3]-t[0])
p[1]=p[0]
fld.update(CP.PT_INPUTS, p[1], t[1])
h[1]=fld.hmass()

"punto 4"
h[4]=h[3]-h[1]+h[0]
p[4]=p[3]
fld.update(CP.HmassP_INPUTS, h[4], p[4])
t[4]=fld.T()

"punto 5"
p[5]=p[0]
h[5]=h[4]
fld.update(CP.HmassP_INPUTS, h[5], p[5])
t[5]=fld.T()
       
"punto 2"
t[2]=c.Temperatura_mandata(t[0]-T_ref, t[1]-T_ref, p[3]*10**-5)+T_ref
p[2]=p[3]
fld.update(CP.PT_INPUTS, p[2], t[2])
h[2]=fld.hmass()


cop=(h[1]-h[3])/(h[2]-h[1])
print('cop senza olio='+str(cop))


plt.figure(dpi=200)
gt.grafico_PH_semplice(P/100000,H/1000, 'r',0)
gt.grafico_PH_semplice(p/100000,h/1000, 'b',1)
plt.plot(H[0]/1000,P[0]/100000,'r',label='con olio\ncop= '+str(round(COP, 2)))    
plt.plot(H[0]/1000,P[0]/100000,'b',label='senza olio\ncop= '+str(round(cop, 2)))    

plt.legend(loc='lower right')
plt.title('Diagramma P-H Thome')