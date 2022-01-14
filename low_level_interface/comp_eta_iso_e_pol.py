"""
eta_pol
"""
import numpy as np
import compressore as c
from CoolProp.CoolProp import PropsSI

lib="REFPROP::"
fluid="CO2"

T=np.zeros(3)
P=np.zeros(3)
H=np.zeros(3)
S=np.zeros(3)
V=np.zeros(3)


T_ref=273.15

T[0]=-5+T_ref
P[0]=PropsSI('P','T',T[0],'Q',1,lib+fluid)

T[1]=30+T_ref
P[1]=P[0]
H[1]=PropsSI('H','T',T[1],'P',P[1],lib+fluid)
S[1]=PropsSI('S','T',T[1],'P',P[1],lib+fluid)
V[1]=1/PropsSI('D','T',T[1],'P',P[1],lib+'CO2')


P[2]=95*10**5
T[2]=c.Temperatura_mandata(T[0]-T_ref, T[1]-T_ref, P[2]*10**-5)+T_ref
H[2]=PropsSI('H','T',T[2],'P',P[2],lib+fluid)
V[2]=1/PropsSI('D','T',T[2],'P',P[2],lib+'CO2')
H2_id=PropsSI('H','S',S[1],'P',P[2],lib+fluid)


R=PropsSI(fluid,"gas_constant") #molare
R=R/PropsSI(fluid,"M")


# =============================================================================
# R1=PropsSI("Cpmass","T",T[1],"P", P[1],lib+"CO2")-PropsSI("Cvmass","T",T[1],"P", P[1],lib+"CO2")
# R2=PropsSI("Cpmass","T",T[2],"P", P[2],lib+"CO2")-PropsSI("Cvmass","T",T[2],"P", P[2],lib+"CO2")
# R=R2
# 
# R=P[1]/PropsSI('D','T',T[1],'P',P[1],lib+'CO2')/T[1]
# =============================================================================

beta=P[2]/P[1]

def pol(n):
    L=n/(n-1)*R*T[1]*(beta**((n-1)/n)-1)
    return L*PropsSI('Z','T',T[1],'P',P[1],lib+fluid)

def f_pol(n): #1.548
    L=pol(n)
    return H[2]-H[1]-L

def f_iso(n):#1.275
    L=pol(n)
    return H2_id-H[1]-L

gamma1=PropsSI("isentropic_expansion_coefficient","T",T[1],"P", P[1],lib+fluid)
gamma2=PropsSI("Cpmass","T",T[1],"P", P[1],lib+fluid)/PropsSI("Cvmass","T",T[1],"P", P[1],lib+fluid)
gamma3=PropsSI("Cp0mass","T",T[1],"P", P[1],lib+fluid)/PropsSI("Cvmass","T",T[1],"P", P[1],lib+fluid)




"PV^n=cost"


n=np.log(beta)/np.log(V[1]/V[2]) #1.384659932187524  #1.385





