import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
import compressore as c
#from mixture_impianto_senza_eiettore_sep_function_T7fix import Funz as Funz7
from mixture_impianto_senza_eiettore_sep_function import Funz
import matplotlib.pyplot as plt
from joblib import Parallel, delayed



lib="REFPROP"

fluid1="CO2&"
fluid2=[] ; x0=[]

fluid2.append("R1234YF") ; x0.append(0.9112132132132) 
fluid2.append("R1233ZD") ; x0.append(0.98) 
#fluid2.append("R1233ZD") ; x0.append(0.985) #se T7=-20


fluid2.append("R1234ZEE"); x0.append(0.98) 
fluid2.append("PENTANE" ); x0.append(0.9955) 
#fluid2.append("HEXANE" ) ; x0.append(0.99) 
fluid2.append("PROPANE" ) ; x0.append(0.9) 
fluid2.append("ISOBUTANE" ) ; x0.append(0.97) 
#fluid2.append("dimethylether" ) ; x0.append(0.97) 
#fluid2.append("CO2"); x0.append(0.9)  #cop=2.25

fluids_list=[fluid1+e for e in fluid2]



n=20
eps=0.8
T_gc=40
#10**np.linspace(0,np.log10(100),10)    =   np.logspace(0,np.log10(100),10)
P_gc=np.linspace(70,110,n)#95
#P_gc=np.linspace(110,70,n)#95
cop=np.zeros(n)
x=np.zeros(n)
T_eva=-5
T_sep=10
T_ref=273.15
eta_c=0.8
xx=[]
copcop=[]
P_gcP_gc=[]


glide=15#10
" itero sulla x finchè T_7 non rispetta il glide"



xx=np.load('xx_opt.npy')
copcop=np.load('copcop_opt.npy')

# =============================================================================
# 
# 
# fig, (ax0, ax1) = plt.subplots(1, 2,dpi=200, figsize=(8,3))
# for j,fluids in enumerate(fluids_list):
#     ax0.plot(P_gc,copcop[j])
#     ax1.plot(P_gc,np.array(xx[j])*100,label=fluid2[j])
# #ax0.plot([90,90],[0,0])    
# ax1.set_xlabel("$ P_{gc} $ [bar]")
# ax0.set_xlabel("$ P_{gc} $ [bar]")
# 
# ax0.set_ylabel("COP []")
# ax1.set_ylabel("x [%$kg_{CO_{2}}/kg_{tot}$]")
# 
# #ax0.set_title('Curva caratteristica $ P_{gc} $, $ T_{gc} $=40°C, $ T_{eva} $='+str(T_eva-glide)+'°C, glide='+str(glide)+'°C')
# ax0.grid()
# ax1.grid()
# #ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.9),ncol=3) 
# #ax1.legend(loc='lower center', bbox_to_anchor=(-0.1, -0.3),ncol=3)    
# ax1.legend(loc='lower center', bbox_to_anchor=(-0.25, -0.6),ncol=3)    
# 
#    
# fig.suptitle('Curva caratteristica $ P_{gc} $, $ T_{gc} $=40°C, $ T_{eva} $='+str(T_eva-glide)+'°C, glide='+str(glide)+'°C')
# #plt.tight_layout()
# plt.subplots_adjust(#left=0.1,
#                     #bottom=0.1, 
#                     #right=0.9, 
#                     #top=0.9, 
#                     wspace=0.35, 
#                     #hspace=0.4
#                     )
# 
# =============================================================================



"ottimi"
opt_index=[]
P_opt=[]
x_opt=[]
for j,fluids in enumerate(fluids_list):
    index=np.argmax(copcop[j])
    opt_index.append(index)
    P_opt.append(P_gc[index])
    x_opt.append(xx[j,index])
 
    
P8=[]
    
plt.figure(dpi=200)   
#plt.plot((0,1),(T_eva-glide+5-T_ref,T_eva+5-T_ref),'r') 
for j,fluids in enumerate(fluids_list):
    print('\n'*10+'fluids = '+fluids)

    mix   = CP.AbstractState(lib, fluids)
    mix_l = CP.AbstractState(lib, fluids)
    mix_g = CP.AbstractState(lib, fluids)

    mix.set_mass_fractions([x_opt[j], 1-x_opt[j]])
  
    funz=Funz(eps,P_opt[j],T_gc,T_eva,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    
    print(cop)
    
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    T,P,H,S=funz.entropie(T,P,H,S)
    
    P8.append(P[8])

    "evaporatore"
    n=100
    H_e=np.linspace(H[7],H[8],n)
    T_e=np.zeros(n)
    for i in range(n):
        mix_l.update(CP.HmassP_INPUTS, H_e[i], P[7])
        T_e[i]=mix_l.T()

    plt.plot(H_e,T_e-T_ref,label=fluid2[j])


plt.xlabel("H [kJ/kg]")
plt.ylabel("T [°C]")
plt.title('Evaporatore')
plt.legend()
plt.grid()





plt.figure(dpi=200)   
#plt.plot((0,1),(T_eva-glide+5-T_ref,T_eva+5-T_ref),'r') 
for j,fluids in enumerate(fluids_list):
    print('\n'*10+'fluids = '+fluids)

    mix   = CP.AbstractState(lib, fluids)
    mix_l = CP.AbstractState(lib, fluids)
    mix_g = CP.AbstractState(lib, fluids)

    mix.set_mass_fractions([x_opt[j], 1-x_opt[j]])
  
    funz=Funz(eps,P_opt[j],T_gc,T_eva,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    
    print(cop)
    
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    T,P,H,S=funz.entropie(T,P,H,S)

    "evaporatore"
    n=100
    H_e=np.linspace(H[7],H[8],n)
    T_e=np.zeros(n)
    for i in range(n):
        mix_l.update(CP.HmassP_INPUTS, H_e[i], P[7])
        T_e[i]=mix_l.T()

    plt.plot(np.linspace(0,100,n),T_e-T_ref,label=fluid2[j])


plt.xlabel("Q [%]")
plt.ylabel("T [°C]")
plt.title('Evaporatore')
plt.legend()

plt.grid()