import numpy as np 
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import compressore as c
import grafici_termodinamici as gt



# =============================================================================
# "Predrag Hrnjak"
# a1 = 6.59664
# a2 = -1.61705E+03
# a3 = 7.00550E+04
# a4 = 3.88634
# a5 = -1.46846E+03
# a6 = 1.19438E+05
# a7 = 3.40692E-01
# a8 = -2.54276E+02
# a9 =  2.11410E+04
# 
# plt.figure(dpi=200)
# plt.grid()
# T=np.linspace(0,10)+273.15
# omega=0.05#np.linspace(0.05,0.4,5)
# #np.log10(P) = a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2
# P = 10**(a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2)
# plt.plot(T,P)
# 
# omega=0.1#np.linspace(0.05,0.4,5)
# #np.log10(P) = a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2
# P = 10**(a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2)
# plt.plot(T,P)
# 
# omega=0.2#np.linspace(0.05,0.4,5)
# #np.log10(P) = a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2
# P = 10**(a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2)
# plt.plot(T,P)
# 
# omega=0.4#np.linspace(0.05,0.4,5)
# #np.log10(P) = a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2
# P = 10**(a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2)
# plt.plot(T,P)
# 
# omega=0.99#np.linspace(0.05,0.4,5)
# #np.log10(P) = a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2
# P = 10**(a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2)
# plt.plot(T,P)
# 
# omega=1#np.linspace(0.05,0.4,5) #no il 100% lo fa praticamente al 40
# #np.log10(P) = a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2
# P = 10**(a1 + a2/T + a3/(T**2) + (a4 + a5/T + a6/(T**2))*np.log10(omega) + (a7 + a8/T + a9/(T**2))*np.log10(omega)**2)
# plt.plot(T,P)
# =============================================================================



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




def P_s(w,T_bub): #bisogna essere dentro la campana

    P_c = (a1+a2*w + a3*w**2)/(1 + a4*w + a5*w**2)  #pressione pseudo critica in bar
    T_c = (b1 + b2*w)/(1 + b3*w + b4*w**2)
    alpha = (c1+c2*w)/(1+c3*w)


    tau=np.log(T_c/T_bub)

    #np.log(P_c/P_s)= alpha*tau + 5.957*tau**2.64
    
    #np.log(P_s) = np.log(P_c) - alpha*tau - 5.957*tau**2.64
    P_ss=np.exp( np.log(P_c) - alpha*tau - 5.957*tau**2.64 )
    
    return P_ss*10**5

def f_P_s(T_bub,w):
    #print(T_bub)
    #print(w)
    #print(P_s(w,T_bub))
    return P_s(w,T_bub)-P[0]
    
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

def P_3(p3):
    P[3]=p3
    
    oil.update(CP.PT_INPUTS, P[3], T[3])
    H_o[3]=oil.hmass()

    H_f[3]=( H[3]- k*H_o[3] )/(1-k)
    
    fld.update(CP.HmassT_INPUTS, H_f[3], T[3])
    Q[3]=fld.Q()

# =============================================================================
#     fld.update(CP.QT_INPUTS, 0, T[3])
#     H_f_l = fld.hmass()
# 
# 
#     fld.update(CP.QT_INPUTS, 1, T[3])
#     H_f_v = fld.hmass()
# 
#     Q[3] = ( H_f[3] - H_f_l )/(H_f_v-H_f_l)
# =============================================================================

# =============================================================================
# fld.update(CP.HmassT_INPUTS, H_f[2], T[3])
# fld.Q()
# fld.p()
# =============================================================================

    z_l=(1-k-Q[3]+k*Q[3])/(1-Q[3]+k*Q[3])

    w[3]=1-z_l

    return P_s(w[3],T[3])-p3

P[3]=fsolve(P_3, P[3])


"punto 0"
P[0]=P[3]
w[0]= 0.56  #è ottimizzabile siamo noi a deciderlo
Q_tot=1-k/w[0]
z_l=1-w[0]
Q[0]=(1-z_l-k)/(1-z_l-k+z_l*k)

T[0]=T[3] #inizializzo

T[0]=fsolve(f_P_s,T[0],args=(w[0]))


oil.update(CP.PT_INPUTS, P[0], T[0])
H_o[0]=oil.hmass()

# =============================================================================
# fld.update(CP.QT_INPUTS, 0, T[0])
# H_f_l = fld.hmass()
# 
# 
# fld.update(CP.QT_INPUTS, 1, T[0])
# H_f_v = fld.hmass()
# 
# H_f[0] = H_f_l*(1-Q[0]) + H_f_v*Q[0]
# =============================================================================

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


COP=(H[0]-H[2])/(H[1]-H[0])/Q[0]
print('cop con olio='+str(COP))


# =============================================================================
# P=P/100000
# H=H/1000
# 
# plt.figure(dpi=200)
# plt.grid()
# plt.xlabel('H [kJ/kg]')
# plt.ylabel('P [bar]')
# plt.plot(H,P,'r')
# plt.plot([H[3],H[0]],[P[3],P[0]],'r')
# =============================================================================







# =============================================================================
# plt.figure(dpi=200)
# plt.grid()
# plt.xlabel('T [°C]')
# plt.ylabel('P [bar]')
# 
# n=6
# m=50
# ww=np.linspace(0,0.6,n)
# TT=np.linspace(-20,30,m)+T_ref
# PP=np.zeros(m)
# 
# for i in range(n):
#     plt.plot(TT-T_ref,P_s(ww[i],TT)/100000)
# 
# for i in range(m):
#     fld.update(CP.QT_INPUTS, 0, TT[i])
#     PP[i]=fld.p()
# 
# plt.plot(TT-T_ref,PP/100000,'--')
# 
# fld.p_critical()
# fld.T_critical()-273.15
# =============================================================================



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