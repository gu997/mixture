import numpy as np 
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP


a=[-2394.5, 182.52, -724.21, 3868.0, -5268.9]
b=[8.0736, -0.72212, 2.3914, -13.779, 17.066]

#w=0
#P_sat=0.55#la vuole in MPa

def T_bub(w,P_sat):
    A=a[0]+a[1]*w+a[2]*w**3+a[3]*w**5+a[4]*w**7
    B=b[0]+b[1]*w+b[2]*w**3+b[3]*w**5+b[4]*w**7

    return A/(np.log(P_sat)-B)

# =============================================================================
# n=100
# w=np.linspace(0,0.6,n)
# T=np.zeros(n)
# 
# T=T_bub(w,P_sat)
# 
# plt.figure(dpi=200)
# plt.plot(w,T-273.15)
# =============================================================================

lib="REFPROP"

fld   = CP.AbstractState(lib, "R134A")

P_sat=3.431*10**5
P1=P_sat+100
P2= P_sat-100
fld.update(CP.PQ_INPUTS, P1 , 0)
#fld.p()
#fld.Q()
T1= fld.T()
#fld.hmass()

fld.update(CP.PQ_INPUTS, P2 , 0)
#fld.p()
#fld.Q()
T2= fld.T()
#fld.hmass()

"trovo a0 e b0 considerando w=0"
P1=P1*10**-6
P2=P2*10**-6

b[0]=(T2*np.log(P2)-T1*np.log(P1))/(T2-T1)

a[0]=T1*(np.log(P1)-b[0])

# =============================================================================
# w=0
# T_bub(w,P_sat*10**-6)-273.15
# =============================================================================

w=0.03

"adesso divido l'evaporazione in 11 punti"
n=11
Q_in=0.15
Q_out=0.95
Q=np.linspace(Q_in,Q_out,n)
T_b=np.zeros(n)
dh=np.zeros(n-1)
h_l=np.zeros(n)
cp_l_ref=np.zeros(n)
cp_l=np.zeros(n)
h_v=np.zeros(n)
cp_v=np.zeros(n)
h=np.zeros(n)




w_loc=w/(1-Q)
T_b=T_bub(w_loc,P_sat*10**-6)

"calcolo calore latente e cp"
for i in range(n):
    fld.update(CP.QT_INPUTS, 0, T_b[i])
    h_l[i]=fld.hmass()
    cp_l_ref[i]=fld.cpmass()
    
    fld.update(CP.QT_INPUTS, 1, T_b[i])
    h_v[i]=fld.hmass()
    cp_v[i]=fld.cpmass()
    
h_lv=h_v-h_l

"calcolo cp oil"
# =============================================================================
# oil   = CP.AbstractState("INCOMP", "SAB")
# oil.update(CP.PQ_INPUTS, P_sat , T_b[i])
# =============================================================================
rho_oil=971
s=rho_oil/999
cp_oil=4.186*(0.388+0.00045*(1.8*T_b+32))/(s**1/2)

cp_l=w_loc*cp_oil+(1+w_loc)*cp_l_ref

for i in range(n-1):
    dh[i]=h_lv[i]*(Q[i+1]-Q[i])+(1-Q[i])*(T_b[i+1]-T_b[i])*cp_l[i]+(T_b[i+1]-T_b[i])*cp_v[i]
    
    h[i+1]=h[i]+dh[i]

dh_tot=h[-1]

"senza olio"
fld.update(CP.PQ_INPUTS, P_sat, 0)
h0=fld.hmass()
fld.update(CP.PQ_INPUTS, P_sat, 1)
h1=fld.hmass()
dh_s_o=h1-h0

