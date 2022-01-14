import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c
from mixture_impianto_senza_eiettore_sep_function import Funz



lib="REFPROP"

fluids="CO2&R1234YF"

mix   = CP.AbstractState(lib, fluids)
mix_l = CP.AbstractState(lib, fluids)
mix_g = CP.AbstractState(lib, fluids)

x=0.6
#mix.set_mole_fractions([x, 1-x])
mix.set_mass_fractions([x, 1-x])


eps=0.8
T_gc=40
P_gc=95
T_eva=-5
T_sep=10
T_ref=273.15



funz=Funz(eps,P_gc,T_gc,T_eva,T_sep,mix,mix_l,mix_g)

T,P,H,S,m,cop=funz.imp()
print(cop)





# =============================================================================
# T,P,H,S=funz.punti_rimanenti(T,P,H,S)
# gt.grafico_PH_sep_IHX(P/100000,H/1000, mix, mix_l, mix_g)
# =============================================================================
