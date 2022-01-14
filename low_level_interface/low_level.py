#import CoolProp
import CoolProp.CoolProp as CP
mix   = CP.AbstractState("REFPROP", "R1234YF&R32&CO2") #HEOS
mix_l = CP.AbstractState("REFPROP", "R1234YF&R32&CO2")
mix_g = CP.AbstractState("REFPROP", "R1234YF&R32&CO2")


mix.set_mole_fractions([0.57897045, 0.36141605, 0.05961349])
P_tot=30*10**5
Q=0.6
mix.update(CP.PQ_INPUTS, P_tot, Q)
#mix.update(CoolProp.PT_INPUTS, 517685.170540156, 278.15)
#mix.update(CP.QT_INPUTS, 0.27193457276482425, 278.15)
m_l=mix.mole_fractions_liquid()
m_g=mix.mole_fractions_vapor()


mix.p()
mix.Q()
T=mix.T()
mix.hmass()


mix_l.set_mole_fractions(m_l)
mix_l.update(CP.PT_INPUTS, P_tot, T)

mix_l.p()
mix_l.Q()
mix_l.T()
mix_l.hmass()

mix_g.set_mole_fractions(m_g)
mix_g.update(CP.PT_INPUTS, P_tot, T)

mix_g.p()
mix_g.Q()
mix_g.T()
mix_g.hmass()


