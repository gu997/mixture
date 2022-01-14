import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
#from mixture_impianto_senza_eiettore_sep_function_T7fix import Funz as Funz7
from mixture_impianto_doppio_evap_function_T6 import Funz
import matplotlib.pyplot as plt
from joblib import Parallel, delayed



lib="REFPROP"

fluid1="CO2&"
fluid2=[] ; x0=[]

fluid2.append("R1234YF") ; x0.append(0.90) 
# =============================================================================
# fluid2.append("R1233ZD") ; x0.append(0.98) 
# #fluid2.append("R1233ZD") ; x0.append(0.985) #se T7=-20
# 
# 
# fluid2.append("R1234ZEE"); x0.append(0.98) 
# fluid2.append("PENTANE" ); x0.append(0.9955) 
# #fluid2.append("HEXANE" ) ; x0.append(0.998) 
# fluid2.append("PROPANE" ) ; x0.append(0.9) 
# fluid2.append("ISOBUTANE" ) ; x0.append(0.97) 
# fluid2.append("dimethylether" ) ; x0.append(0.97) 
# =============================================================================
fluid2.append("CO2"); x0.append(0.99)  #cop=2.25

fluids_list=[fluid1+e for e in fluid2]



n=20
eps=0.8
T_gc=40
#10**np.linspace(0,np.log10(100),10)    =   np.logspace(0,np.log10(100),10)
#P_gc=np.linspace(70,110,n)#95
P_gc=np.linspace(110,70,n)#95
cop=np.zeros(n)
x=np.zeros(n)
T_acqua=-5
#T_EVA_AT=-5
#T_eva=-5
T_ref=273.15
eta_c=0.8
eta1=0.9
eta2=0.8
eta_mix=0.85
v=10 #m/s
xx=[]
copcop=[]
P_gcP_gc=[]

dT=10
" itero sulla x finchè T_7 non rispetta il glide"



def f_T_EVA_AT_CO2(y,P_gcc):
    f=np.zeros(1)
    
    T_eva=T_acqua-dT
    T_EVA_AT=y
    x=1
    mix.set_mass_fractions([x, 1-x])
    
    funz=Funz(eps,P_gcc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
    
    T,P,H,S,m,cop=funz.imp()
    
    xx=(H[7]-H[6])*m[2]/(H[0]-H[9]+(H[7]-H[6])*m[2])
    
    #f[0]=(T[9]-(dT*xx+T_eva+T_ref))
    #f[1]=(T[7]-(dT*xx+T_eva+T_ref))
    
    f[0]=(T[9]-(dT*xx+T_acqua-dT+T_ref))
    #f[1]=(T[7]-(dT*xx+T_acqua-dT+T_ref))
    
    print('constr= ',f)
    print(y)
    return f

def f_T_EVA_AT(y,P_gcc):
    f=np.zeros(2)
    
    T_eva=y[0]
    T_EVA_AT=T_acqua
    x=y[1]
    mix.set_mass_fractions([x, 1-x])
    
    funz=Funz(eps,P_gcc,T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
    
    T,P,H,S,m,cop=funz.imp()
    
    xy=(H[7]-H[6])*m[2]/(H[0]-H[9]+(H[7]-H[6])*m[2])
    
    #f[0]=(T[9]-(dT*xx+T_eva+T_ref))
    #f[1]=(T[7]-(dT*xx+T_eva+T_ref))
    
    f[0]=(T[9]-(dT*xy+T_EVA_AT-dT+T_ref))
    f[1]=(T[7]-(dT*xy+T_EVA_AT-dT+T_ref))
    
    print('constr= ',f)
    print(y)
    return f
    
    


for j,fluids in enumerate(fluids_list):
    print('\n'*10+'fluids = '+fluids)

    mix   = CP.AbstractState(lib, fluids)
    mix_l = CP.AbstractState(lib, fluids)
    mix_g = CP.AbstractState(lib, fluids)
    
    
    #mix.set_mass_fractions([x, 1-x])
    
    x00=x0[j]
    
    for i in range(n):
        print('i =',i)
        try:
            if fluid2[j]=='CO2':
                y=[-5]
                y=fsolve(f_T_EVA_AT_CO2,y,args=(P_gc[i]))        
                T_EVA_AT=y
                T_eva=T_acqua-dT
                x[i]=1
                mix.set_mass_fractions([x[i], 1-x[i]])
                funz=Funz(eps,P_gc[i],T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
                    
                T,P,H,S,m,cop[i]=funz.imp()
                
                
            
            else:
                y=[-15,x00]
                y=fsolve(f_T_EVA_AT,y,args=(P_gc[i]))        
                T_eva=y[0]
                x[i]=y[1] 
                mix.set_mass_fractions([x[i], 1-x[i]])
                
                T_EVA_AT=T_acqua
                funz=Funz(eps,P_gc[i],T_gc,T_eva,T_EVA_AT,eta_c,eta1,eta2,eta_mix,v,mix,mix_l,mix_g)
                    
                T,P,H,S,m,cop[i]=funz.imp()
            
                #x00=x[i]
                
        #... run your code
        #except:
         #   pass
        except Exception as e:
            print(e)
            x[i]=1
            cop[i]=0
       
        

    print(x)
    
    c=np.where(cop>0)
    copa=cop[c]
    xa=x[c]
    p_gca=P_gc[c]
    
    xl=xa.tolist()
    copl=copa.tolist()
    p_gcl=p_gca.tolist()

    xx.append(xl)
    copcop.append(copl)
    P_gcP_gc.append(p_gcl)

# =============================================================================
#     ax0.plot(P_gc,cop)
#     ax1.plot(P_gc,x,label=fluid2[j])
# =============================================================================
    #ax1.plot(P_gc,T_8,label=fluid2[j])

xx=np.array(xx)
copcop=np.array(copcop)
P_gcP_gc=np.array(P_gcP_gc)

# =============================================================================
# xx=np.load('xx_opt.npy')
# copcop=np.load('copcop_opt.npy')
# =============================================================================


fig, (ax0, ax1) = plt.subplots(1, 2,dpi=200, figsize=(8,3))
for j,fluids in enumerate(fluids_list):
    ax0.plot(P_gcP_gc[j],copcop[j],'--')
    ax1.plot(P_gcP_gc[j],np.array(xx[j])*100,label=fluid2[j])
#ax0.plot([90,90],[0,0])    
ax1.set_xlabel("$ P_{gc} $ [bar]")
ax0.set_xlabel("$ P_{gc} $ [bar]")

ax0.set_ylabel("COP []")
ax1.set_ylabel("x [%$kg_{CO_{2}}/kg_{tot}$]")

#ax0.set_title('Curva caratteristica $ P_{gc} $, $ T_{gc} $=40°C, $ T_{eva} $='+str(T_eva-glide)+'°C, glide='+str(glide)+'°C')
ax0.grid()
ax1.grid()
#ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.9),ncol=3) 
#ax1.legend(loc='lower center', bbox_to_anchor=(-0.1, -0.3),ncol=3)    
ax1.legend(loc='lower center', bbox_to_anchor=(-0.25, -0.6),ncol=3)    

   
fig.suptitle('Curva caratteristica $ P_{gc} $, $ T_{gc} $=40°C, $ T_{acqua_{in}} $='+str(T_acqua)+'°C, glide='+str(dT)+'°C')
#plt.tight_layout()
plt.subplots_adjust(#left=0.1,
                    #bottom=0.1, 
                    #right=0.9, 
                    #top=0.9, 
                    wspace=0.35, 
                    #hspace=0.4
                    )