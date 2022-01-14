#import CoolProp
import CoolProp.CoolProp as CP
#mix   = CP.AbstractState("REFPROP", "R1234YF&R32&CO2")
mix   = CP.AbstractState("REFPROP", "CO2")
mix   = CP.AbstractState("REFPROP", "POE")
CP.set_reference_state('REFPROP::POE','IIR')

mix   = CP.AbstractState("REFPROP", "PENTANE - Copia")

#mix   = CP.AbstractState("HEOS", "T66")


#mix.set_mole_fractions([0.57897045, 0.36141605, 0.05961349])
mix   = CP.AbstractState("REFPROP", "POE&CO2")
CP.set_reference_state('REFPROP::POE','IIR')

mix.set_mass_fractions([0.5 , 0.5])


P1=30*10**5
P2=90*10**5

Q=0.6
T1=273.15
T2=T1+100
#mix.update(CP.PQ_INPUTS, P_tot, Q)
mix.update(CP.PT_INPUTS, P1, T1)
print('mix.hmass()=',mix.hmass())
print('mix.cpmass()',mix.cpmass())

mix.update(CP.PT_INPUTS, P2, T2)
print('mix.hmass()=',mix.hmass())
print('mix.cpmass()',mix.cpmass())

#mix.update(CP.QT_INPUTS, 0.27193457276482425, 278.15)


# =============================================================================
# mix.p()
# mix.Q()
# mix.T()
# mix.hmass()
# mix.cpmass()
# 
# from CoolProp.CoolProp import PropsSI
# PropsSI('H','T',300,'P',101325,'SAB')
# PropsSI('H','T',300,'P',101325,'INCOMP::SAB')
# mix   = CP.AbstractState("INCOMP", "SAB")
# mix.update(CP.PT_INPUTS, 101325, 300)
# mix.p()
# mix.Q()
# mix.T()
# mix.hmass()
# 
# 
# a1=0.10288132e+02    
# a2=-0.26953770e-01    
# a3= 0.20951065e-03    
# a4=-0.27910773e-06
# a5= 0.12266269e-09    
# =============================================================================

#a1+a2*T+a3*T**2+a4*T**3+a5*T**4  +8.3143


# =============================================================================
# jj=CP.get_fluid_param_string("CO2", "JSON")
# 
# file_uno = open("CO2.txt", "w")
# file_uno.write(jj)
# file_uno.close()   
# 
# =============================================================================
a1= 30.0 
a2=0.0
a3=69.0 
a4=850.0
a5=98.0 
a6=2000.0