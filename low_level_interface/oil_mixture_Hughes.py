import CoolProp.CoolProp as CP
import numpy as np 

lib="REFPROP"
fld   = CP.AbstractState(lib, "R12")
oil   = CP.AbstractState("INCOMP", "SAB")

T=np.zeros(4)
P=np.zeros(4)
H=np.zeros(4)

T_ref=273.15


"olio percentuale nella fase liquida"
w=0.08

"punto 2"
T[2]=50+T_ref
fld.update(CP.QT_INPUTS, 0, T[2])
P[2]=fld.p()

oil.update(CP.PT_INPUTS, P[2] , T[2])

H[2] = w*oil.hmass() + (1-w)*fld.hmass()


"punto 3"
T[3]=0+T_ref
H[3]=H[2]

fld.update(CP.QT_INPUTS, 1, T[3])
h_vap=fld.hmass()

fld.update(CP.QT_INPUTS, 0, T[3])
h_r_liq=fld.hmass()
pp=fld.p()

oil.update(CP.PT_INPUTS, pp , T[3])
h_oil=oil.hmass()

h_liq=w*h_oil+(1-w)*h_r_liq

#h_mix=Q*h_vap + (1-Q)*h_liq
Q=(H[3]-h_liq)/(h_vap-h_liq)



