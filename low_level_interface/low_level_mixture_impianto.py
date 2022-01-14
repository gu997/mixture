import numpy as np
import CoolProp.CoolProp as CP
import grafici_termodinamici_mixture as gt

lib="REFPROP"

mix   = CP.AbstractState(lib, "R1234YF&R32&CO2")
mix_l = CP.AbstractState(lib, "R1234YF&R32&CO2")
mix_g = CP.AbstractState(lib, "R1234YF&R32&CO2")

#mix.set_mole_fractions([0.57897045, 0.36141605, 0.05961349])
mix.set_mass_fractions([0.755, 0.215, 0.03])



T_gc=40
T_eva=-0.88#5
T_sep=5
eta_c=0.8


T=np.zeros(9)
P=np.zeros(9)
H=np.zeros(9)
S=np.zeros(9)
m=np.ones(3)

T_ref=273.15

"PUNTO 2"

T[2]=T_gc+T_ref
mix.update(CP.QT_INPUTS, 0, T[2])
P[2]=mix.p()
H[2]=mix.hmass()
S[2]=mix.smass()

"PUNTO 3"
T[3]=T_sep+T_ref 
H[3]=H[2]
mix.update(CP.HmassT_INPUTS, H[3], T[3])
P[3]=mix.p()
S[3]=mix.smass()
Q=mix.Q()



"SEPARATORE"
m_l=mix.mole_fractions_liquid()
m_g=mix.mole_fractions_vapor()


mix_l.set_mole_fractions(m_l)
mix_g.set_mole_fractions(m_g)
#Q=Q*mix_g.molar_mass()/mix.molar_mass()

"PUNTO 4"
mix_l.update(CP.PT_INPUTS, P[3], T[3])
P[4]=mix_l.p()
mix_l.Q()
T[4]=mix_l.T()
H[4]=mix_l.hmass()
S[4]=mix_l.smass()

"PUNTO 6"
T[6]=T_eva+T_ref
mix_l.update(CP.QT_INPUTS, 1, T[6])
P[6]=mix_l.p()
H[6]=mix_l.hmass()
S[6]=mix_l.smass()

"PUNTO 5"
H[5]=H[4]
P[5]=P[6]
mix_l.update(CP.HmassP_INPUTS, H[5], P[5])
T[5]=mix_l.T()
S[5]=mix_l.smass()

"PUNTO 7"
T[7]=T[3]
P[7]=P[3]
mix_g.update(CP.PT_INPUTS, P[3], T[3])
mix_g.Q()
H[7]=mix_g.hmass()
S[7]=mix_g.smass()

"PUNTO 8"
H[8]=H[7]
P[8]=P[6]
mix_g.update(CP.HmassP_INPUTS, H[8], P[8])
T[8]=mix_g.T()
S[8]=mix_g.smass()

"PUNTO 0"
Q=Q*mix_g.molar_mass()/mix.molar_mass()
m[1]=1-Q
m[2]=Q
H[0]=m[1]*H[6]+m[2]*H[8]
P[0]=P[6]
mix.update(CP.HmassP_INPUTS, H[0], P[0])
T[0]=mix.T()
S[0]=mix.smass()

"PUNTO 1"
P[1]=P[2]
mix.update(CP.PSmass_INPUTS, P[1], S[0])
H1_id=mix.hmass()
H[1]=H[0]+(H1_id-H[0])/eta_c
mix.update(CP.HmassP_INPUTS, H[1], P[1])
T[1]=mix.T()
S[1]=mix.smass()


cop=m[1]*(H[6]-H[5])/(H[1]-H[0])
cop2=(H[0]-H[2])/(H[1]-H[0])
print('cop =',cop)
print('cop2=',cop2)


gt.grafico_PH_sep(P/100000,H/1000, mix, mix_l, mix_g)

# =============================================================================
# H=H*10**-3
# P=P*10**-5
# plt.plot(H[:7],P[:7],'r')
# plt.plot([H[3],H[7],H[8],H[0]],[P[3],P[7],P[8],P[0]],'r')
# #plot.Show()
# =============================================================================


print(T[5]-T_ref)
