   
import time
import numpy as np
import CoolProp.CoolProp as CP

from CoolProp.CoolProp import PropsSI

def property_calc(args): 
    return PropsSI(*args)

if __name__ == "__main__":
    fluid = 'CO2'
    T_triple = PropsSI(fluid, 'Ttriple')
    T_crit = PropsSI(fluid, 'Tcrit')

    params = [['D', 'T', Ti, 'Q', 1, fluid] for Ti in np.linspace(T_triple, T_crit, 10000)]

    start = time.time()
    results = list(map(property_calc, params))
    print(time.time() - start)



class F:  
    def __init__(self,mix):        
        self.mix = mix    
    def property_calcc(self,args): 
        
        self.mix.update(CP.QT_INPUTS, 1, *args)
        H=mix.hmass()
        S=mix.smass()
        return H#PropsSI('D', 'T', *args, 'Q', 1, 'R717')

if __name__ == "__main__":
    lib="REFPROP"

    #fluids="CO2&R1234YF"
    #fluids="CO2&PROPANE"
    #fluids="CO2&ISOBUTANE"
    #fluids="CO2&HEXANE"
    #fluids="CO2&R1233ZD"


    fluids = 'CO2'
    mix   = CP.AbstractState(lib, fluids)
    f=F(mix)
    
    
    #fluid = 'R717'
    T_triple = PropsSI(fluid, 'Ttriple')
    T_crit = PropsSI(fluid, 'Tcrit')

    Ti=np.linspace(T_triple, T_crit, 10)
    paramss = [[Ti] for Ti in np.linspace(T_triple, T_crit, 10000)]
 

    start = time.time()
    results = list(map(f.property_calcc, paramss))
    print(time.time() - start)
    
    
    
    
import time
import numpy as np
from multiprocessing import Pool

from CoolProp.CoolProp import PropsSI

def property_calc(args): 
    return PropsSI(*args)

if __name__ == "__main__":
   # fluid = 'R717'
    #T_triple = PropsSI(fluid, 'Ttriple')
    #T_crit = PropsSI(fluid, 'Tcrit')

    #params = [['D', 'T', Ti, 'Q', 1, fluid] for Ti in np.linspace(T_triple, T_crit, 5000)]

    start = time.time()
    with Pool(processes=12) as pool:
        results = list(pool.map(property_calc, params))    
    print(time.time() - start)
