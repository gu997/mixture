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

a2_oil = 699.4
a3_oil = 3.976
 
 
 
#a2_oil = 1158.8
#a3_oil = 2.3639

def h_oil(T): #preso da Rémi Dickes

    h_r=a2_oil*273.15 + a3_oil/2*273.15**2
    h=a2_oil*T + a3_oil/2*T**2  -h_r
    
    return h/1000 + 200000


def T_oil(h): #preso da Rémi Dickes
    h_r=a2_oil*273.15 + a3_oil/2*273.15**2
    h=(h-200000)*1000 + h_r
    #h=a2_oil*T + a3_oil/2*T**2
    t=(-a2_oil + (a2_oil**2 + 4*a3_oil/2*h)**(1/2))/a3_oil
    
    return t


def P_s(w,T_bub): #bisogna essere dentro la campana

    P_c = (a1+a2*w + a3*w**2)/(1 + a4*w + a5*w**2)  #pressione pseudo critica in bar
    T_c = (b1 + b2*w)/(1 + b3*w + b4*w**2)
    alpha = (c1+c2*w)/(1+c3*w)


    tau=np.log(T_c/T_bub)

    #np.log(P_c/P_s)= alpha*tau + 5.957*tau**2.64
    
    #np.log(P_s) = np.log(P_c) - alpha*tau - 5.957*tau**2.64
    P_ss=np.exp( np.log(P_c) - alpha*tau - 5.957*tau**2.64 )
    
    return P_ss*10**5



eps=0.8
T_gc=40
P_gc=95
T_eva=-5
#T_sep=10
P_sep=45*10**5#4502182.914250896



T=np.zeros(12)
P=np.zeros(12)
H=np.zeros(12)
H_f=np.zeros(12)
H_o=np.zeros(12)
Q=np.zeros(12)
Q_tot=np.zeros(12)
S=np.zeros(12)
w=np.ones(12)
m=np.ones(3)

T_ref=273.15

lib="REFPROP"
fld   = CP.AbstractState(lib, "CO2")  

"PUNTI NOTI A PRIORI"

"punto 3"
k=0.03#0.0197 #fraz di olio nel totale   massa olio/massa totale
T[3]=T_gc+T_ref
P[3]=P_gc*10**5   #iperparametro
fld.update(CP.PT_INPUTS, P[3], T[3])
H_f[3]=fld.hmass()

H_o[3]=h_oil(T[3])

H[3]=k*H_o[3] + (1-k)*H_f[3]

"punto 8"
T[8]=T_eva+T_ref
w[8]= 0.6  #è ottimizzabile siamo noi a deciderlo
#Q_tot=1-k2/w[8]
#z_l=1-w[8]
#Q[8]=(1-z_l-k2)/(1-z_l-k2+z_l*k2)
P[8]=P_s(w[8],T[8])

H_o[8]=h_oil(T[8])

#fld.update(CP.QT_INPUTS, Q[8], T[8])
#H_f[8]=fld.hmass()

#H[8] = k2*H_o[8] + (1-k2)*H_f[8]


"punto 0"
P[0]=P[8]
#T[0]=T[8]
        
"punto 5"
P[5]=P_sep
#P[5]=PropsSI('P','Q',0,'T',T[5],'CO2')

"punto 11"#come punto 0 ma pura CO2
P[11]=P[8]
Q[11]=1
fld.update(CP.PQ_INPUTS, P[11], Q[11])
H[11]=fld.hmass()
T[11]=fld.T()

"punto 9" #ipotizzo sia CO2 pura
P[9]=P[5]
#"punto 9"
#T[9]=T[5]
fld.update(CP.PQ_INPUTS, P[9], 1)
#fld.update(CP.QT_INPUTS, 1, T[9])
H[9]=fld.hmass() 
T[9]=fld.T()

"punto 10"
P[10]=P[8]
#H[10]=H[11]# così sarebbe uguale al caso senza separatore
H[10]=H[9]
fld.update(CP.HmassP_INPUTS, H[10], P[10])
T[10]=fld.T()    

"punto 6"
P[6]=P[5]
#T[6]=T[5]
#H[6]=PropsSI('H','T',T[5],'Q',0,'CO2')

"punto 7"
#H[7]=H[6]
P[7]=P[8]
#T[7]=PropsSI('T','H',H[7],'P',P[7],'CO2')
  
"punto 10"     
#H[10]=H[9]
#P[10]=P[8]
#T[10]=PropsSI('T','H',H[10],'P',P[10],'CO2')

"punto 1"
#T[1]=T[0] + eps*(T[3]-T[0])
P[1]=P[0]
#H[1]=PropsSI('H','T',T[1],'P',P[1],'CO2')

"punto 4"
P[4]=P[3]

"punto 2"
P[2]=P[3]

"punto 1"
#T[1]=T[0] + eps*(T[3]-T[0])
T[1]=T[11] + eps*(T[3]-T[11])
fld.update(CP.PT_INPUTS, P[1], T[1])
H_f[1]=fld.hmass()
H[1]=H_f[1]
#H_o[1]=h_oil(T[1])

#H[1]=k*H_o[1] + (1-k)*H_f[1]

"IPOTIZZO IL PUNTO 0"
def loop(q):
    "punto 0"
    Q[0]=q
    z_l=(1-k-Q[0]+k*Q[0])/(1-Q[0]+k*Q[0])
    w[0]=1-z_l
    Q_tot[0]=1-k/w[0]
    T[0]=T[8] #inizializzo
    
    def f_P_s(T_bub,w):
        return P_s(w,T_bub)-P[0]
    
    T[0]=fsolve(f_P_s,T[0],args=(w[0]))
    
    fld.update(CP.PQ_INPUTS, P[0], Q[0])
    H_f[0]=fld.hmass()
    
    H_o[0]=h_oil(T[0])
    
    H[0] = k*H_o[0] + (1-k)*H_f[0]
    
    
    
    
    "punto 4"
    #H[4]=H[3]-H[1]+H[0]
    H[4]=H[3]-(H_f[1]-H[11])*Q_tot[0]
    fld.update(CP.HmassP_INPUTS, H[4], P[4])
    T[4]=fld.T()    # T di tentativo
    
    def T_4(t4):
        #print(t4)
    
        T[4]=t4
        
        H_o[4]=h_oil(T[4])
    
        H_f[4]=( H[4]- k*H_o[4] )/(1-k)
        
        fld.update(CP.HmassP_INPUTS, H_f[4], P[4])
    
        return fld.T()-t4
    
    T[4]=fsolve(T_4, T[4])
    
    
    "punto 5"
    H[5]=H[4]
    
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
    
    Q_tot[5]=1-k/w[5]
    m[1]=m[0]*(1-Q_tot[5])
    m[2]=m[0]*Q_tot[5]
    
    "punto 6"
    "bilancio entalpie H[5]=H[9]*Q_tot + H[6]*(1-Q_tot)"
    H[6]=(H[5] - H[9]*Q_tot[5])/(1-Q_tot[5])
    #z_l=(1-k-Q[5]+k*Q[5])/(1-Q[5]+k*Q[5])
    w[6]=w[5]#1-z_l
    k2=w[6]
     
    #T[6]=T[5] #separatore
    
    fld.update(CP.PQ_INPUTS, P[6], 0)
    H_f[6]=fld.hmass()
    
    H_o[6]=(H[6] - (1-k2)*H_f[6])/k2
    T[6]=T[5]#T_oil(H_o[6])
    #H_o[6]=h_oil(T[6])
    
    #H[6]=k2*H_o[6] + (1-k2)*H_f[6]
    
    "punto 7"
    H[7]=H[6]
    
    fld.update(CP.HmassP_INPUTS, H[7], P[7])
    T[7]=fld.T()    # T di tentativo
    
    def T_7(t7):
        T[7]=t7
        
        #oil.update(CP.PT_INPUTS, P[3], T[3])
        H_o[7]=h_oil(T[7])#oil.hmass()
    
        H_f[7]=( H[7]- k2*H_o[7] )/(1-k2)
        
        fld.update(CP.HmassP_INPUTS, H_f[7], P[7])
        Q[7]=fld.Q()
    
        z_l=(1-k2-Q[7]+k2*Q[7])/(1-Q[7]+k2*Q[7])
    
        w[7]=1-z_l
    
        return P_s(w[7],T[7])-P[7]
    
    T[7]=fsolve(T_7, T[7])
    
    "punto 8"
    z_l=1-w[8]
    Q[8]=(1-z_l-k2)/(1-z_l-k2+z_l*k2)
    fld.update(CP.PQ_INPUTS, P[8], Q[8])
    H_f[8]=fld.hmass()
    H[8] = k2*H_o[8] + (1-k2)*H_f[8]
    
# =============================================================================
#     "punto 9"
#     T[9]=T[5]
#     #fld.update(CP.PQ_INPUTS, P[9], 1)
#     fld.update(CP.QT_INPUTS, 1, T[9])
#     H[9]=fld.hmass() 
#     
#     "punto 10"
#     H[10]=H[9]
#     #fld.update(CP.HmassP_INPUTS, H[10], P[10])
#     #T[10]=fld.T()    
# =============================================================================
    
    
    "punto 0"
    
            
    H[0]=(m[1]*H[8]+m[2]*H[10])/m[0]
    fld.update(CP.HmassP_INPUTS, H[0], P[0])
    T[0]=fld.T()    # T di tentativo
    
    def T_0(t0):
        T[0]=t0
        
        #oil.update(CP.PT_INPUTS, P[3], T[3])
        H_o[0]=h_oil(T[0])#oil.hmass()
    
        H_f[0]=( H[0]- k*H_o[0] )/(1-k)
        
        fld.update(CP.HmassP_INPUTS, H_f[0], P[0])
        Q[0]=fld.Q()
    
        z_l=(1-k-Q[0]+k*Q[0])/(1-Q[0]+k*Q[0])
    
        w[0]=1-z_l
    
        return P_s(w[0],T[0])-P[0]
    
    T[0]=fsolve(T_0, T[0])
    
    return Q[0]-q
    
Q[0]=fsolve(loop,0.9)

# =============================================================================
# "punto 10"  
# fld.update(CP.HmassP_INPUTS, H[10], P[10])
# T[10]=fld.T()  
# =============================================================================

"punto 2"
T[2]=c.Temperatura_mandata(T[11]-T_ref, T[1]-T_ref, P[2]*10**-5)+T_ref

fld.update(CP.PT_INPUTS, P[2], T[2])
H_f[2]=fld.hmass()

H_o[2]=h_oil(T[2])

H[2]=k*H_o[2] + (1-k)*H_f[2]

# =============================================================================
# cop=(H[1]-H[3])/(H[2]-H[1])
# print('cop='+str(cop))
# cop=m[1]*(H[8]-H[7])/(H[2]-H[1])
# print('cop='+str(cop))
# =============================================================================

fld.update(CP.PQ_INPUTS, P[0], 0)
H_f_liq_0=fld.hmass()#(1-k)*(1-Q[0])*fld.hmass()#H[0]-H[6]*Q_tot
W_comp=(H_f[2]-H_f[1])*Q_tot[0] + (H_f[2]-H_f_liq_0)*(1-Q_tot[0])*(1-w[0]) + k*(H_o[2]-H_o[0])
# =============================================================================
# print('H_f[2] = ',H_f[2])
# print('H_f[1] = ',H_f[1])
# print('H_f_liq_0 = ',H_f_liq_0)
# print('H_o[2] = ',H_o[2])
# print('H_o[0] = ',H_o[0])
# 
# =============================================================================
m=m*(c.Portata(T[11]-T_ref, T[1]-T_ref, P[2]*10**-5)/3600)
W_comp_c=c.PotenzaE_C(T[11]-T_ref, T[1]-T_ref, P[2]*10**-5)*1000/(m[0]/Q_tot[0])
print('W_comp = ',W_comp)
print('W_comp_c = ',W_comp_c)
print('\n'*3)

COP=(H[2]-H[3]-W_comp)/W_comp
print('COP='+str(COP))

COP=m[1]*(H[8]-H[7])/W_comp/m[0]
print('COP='+str(COP))


gt.grafico_PH_sep(P/100000,H/1000)