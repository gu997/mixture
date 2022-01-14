"considero separatore olio"

import numpy as np 
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import compressore as c
import grafici_termodinamici as gt


"Zhelezny (2007)"
a1 = 73.7537 
a2 = -36.4423 
a3 = -37.0377 
a4 = -1.816879 
a5 = 0.973128
b1 = 303.897 
b2 = -303.666 
b3 = -1.22524 
b4 = 0.225542
c1 = 6.79411 
c2 = -5.6879 
c3 = -0.624642

def h_oil(T): #preso da Rémi Dickes
    a2_oil = 699.4
    a3_oil = 3.976
     
     
     
    #a2_oil = 1158.8
    #a3_oil = 2.3639
    h_r=a2_oil*273.15 + a3_oil/2*273.15**2
    h=a2_oil*T + a3_oil/2*T**2  -h_r
    
    return h/1000 + 200000


def P_s(w,T_bub): #bisogna essere dentro la campana

    P_c = (a1+a2*w + a3*w**2)/(1 + a4*w + a5*w**2)  #pressione pseudo critica in bar
    T_c = (b1 + b2*w)/(1 + b3*w + b4*w**2)
    alpha = (c1+c2*w)/(1+c3*w)


    tau=np.log(T_c/T_bub)

    #np.log(P_c/P_s)= alpha*tau + 5.957*tau**2.64
    
    #np.log(P_s) = np.log(P_c) - alpha*tau - 5.957*tau**2.64
    P_ss=np.exp( np.log(P_c) - alpha*tau - 5.957*tau**2.64 )
    
    return P_ss*10**5


    
T=np.zeros(7)
P=np.zeros(7)
H=np.zeros(7)
H_f=np.zeros(7)
H_o=np.zeros(7)
Q=np.zeros(7)
S=np.zeros(7)
w=np.ones(7)

eps=0.8

T_ref=273.15

lib="REFPROP"

#lib="HEOS"
fld   = CP.AbstractState(lib, "CO2")  
#oil   = CP.AbstractState("INCOMP", "SAB")


"punto 3"  
k=0.03#0.027784 #fraz di olio nel totale   massa olio/massa totale
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
#z_l=1-w[0]
#Q[0]=(1-z_l-k)/(1-z_l-k+z_l*k)
Q[0] = Q_tot/(1-k)
P[0]=P_s(w[0],T[0])

#oil.update(CP.PT_INPUTS, P[0], T[0])
H_o[0]=h_oil(T[0])#oil.hmass()

fld.update(CP.PQ_INPUTS, P[0], Q[0])
H_f[0]=fld.hmass()

H[0] = k*H_o[0] + (1-k)*H_f[0]

"punto 6"#come punto 0 ma pura CO2
P[6]=P[0]
Q[6]=1
fld.update(CP.PQ_INPUTS, P[6], Q[6])
H[6]=fld.hmass()
T[6]=fld.T()



"punto 1"#solo CO2
T[1]=T[6] + eps*(T[3]-T[6]) #in realtà T[0] è solo CO2  riguarda formula keys
#T[1]=T[0] + eps*(T[3]-T[0])
P[1]=P[0]
fld.update(CP.PT_INPUTS, P[1], T[1])
H[1]=fld.hmass()
H_f[1]=fld.hmass()

#oil.update(CP.PT_INPUTS, P[1], T[1])
#H_o[1]=h_oil(T[1])#oil.hmass()

#H[1]=k*H_o[1] + (1-k)*H_f[1]


"punto 4"
H[4]=H[3]-(H_f[1]-H[6])*Q_tot
#H[4]=H[3]-H[1]+H[0]
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
    
    fld.update(CP.HmassP_INPUTS, H_f[5], P[5])
    Q[5]=fld.Q()

    z_l=(1-k-Q[5]+k*Q[5])/(1-Q[5]+k*Q[5])

    w[5]=1-z_l

    return P_s(w[5],T[5])-P[5]

T[5]=fsolve(T_5, T[5])


"punto 2"
P[2]=P[3]
T[2]=c.Temperatura_mandata(T[6]-T_ref, T[1]-T_ref, P[2]*10**-5)+T_ref #nupva T[0]    calcola tutte le H con la PQ

fld.update(CP.PT_INPUTS, P[2], T[2])
H_f[2]=fld.hmass()

#oil.update(CP.PT_INPUTS, P[2], T[2])
H_o[2]=h_oil(T[2])#oil.hmass()

H[2]=k*H_o[2] + (1-k)*H_f[2]


#W_comp=(H_f[2]-H_f[1])*Q[0] + (H_f[2]-H_f[0])*(1-Q[0]) + k*(H_o[2]-H_o[0])
fld.update(CP.PQ_INPUTS, P[6], 0)
H_f_liq_0=fld.hmass()#(1-k)*(1-Q[0])*fld.hmass()#H[0]-H[6]*Q_tot
W_comp=(H_f[2]-H_f[1])*Q_tot + (H_f[2]-H_f_liq_0)*(1-Q_tot)*(1-w[0]) + k*(H_o[2]-H_o[0])
# =============================================================================
# print('H_f[2] = ',H_f[2])
# print('H_f[1] = ',H_f[1])
# print('H_f_liq_0 = ',H_f_liq_0)
# print('H_o[2] = ',H_o[2])
# print('H_o[0] = ',H_o[0])
# =============================================================================

eta=1
m=(c.Portata(T[6]-T_ref, T[1]-T_ref, P[2]*10**-5)/3600)
W_comp_c=c.PotenzaE_C(T[6]-T_ref, T[1]-T_ref, P[2]*10**-5)*1000/(m/Q_tot)*eta
print('W_comp = ',W_comp)
print('W_comp_c = ',W_comp_c)
print('\n'*3)

COP=(H[2]-H[3]-W_comp)/W_comp
#COP=(H[1]-H[3])/(H[2]-H[1])
print('cop con olio='+str(COP))

COP=(H[0]-H[5])/W_comp
#COP=(H[0]-H[5])/(H[2]-H[1])
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

W_comp_co2=c.PotenzaE_C(t[0]-T_ref, t[1]-T_ref, p[3]*10**-5)/(c.Portata(t[0]-T_ref, t[1]-T_ref, p[3]*10**-5)/3600)*1000


cop=(h[1]-h[3])/(h[2]-h[1])
print('cop senza olio=',cop)


plt.figure(dpi=200)
gt.grafico_PH_semplice(P/100000,H/1000, 'r',0)
gt.grafico_PH_semplice(p/100000,h/1000, 'b',1)
plt.plot(H[0]/1000,P[0]/100000,'r',label='con olio\ncop= '+str(round(COP, 2)))    
plt.plot(H[0]/1000,P[0]/100000,'b',label='senza olio\ncop= '+str(round(cop, 2)))    

plt.legend(loc='lower right')
plt.title('Diagramma P-H Zhelezny')
