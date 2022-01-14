import numpy as np
import CoolProp.CoolProp as CP
#import grafici_termodinamici as gt
#import grafici_termodinamici_mixture as gt
from scipy.optimize import fsolve
from mixture_impianto_senza_eiettore_sep_function_T7fix import Funz as Funz7
from mixture_impianto_senza_eiettore_sep_function import Funz
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#from joblib import Parallel, delayed



lib="REFPROP"

fluid1="CO2&"
fluid2=[] ; x0=[]

fluid2.append("R1234YF") ; x0.append(0.9112132132132) 
# =============================================================================
# fluid2.append("R1233ZD") ; x0.append(0.98) 
# #fluid2.append("R1233ZD") ; x0.append(0.985) #se T7=-20
# 
# 
# fluid2.append("R1234ZEE"); x0.append(0.98) 
# fluid2.append("PENTANE" ); x0.append(0.9955) 
# #fluid2.append("HEXANE" ) ; x0.append(0.99) 
# fluid2.append("PROPANE" ) ; x0.append(0.9) 
# fluid2.append("ISOBUTANE" ) ; x0.append(0.97) 
# fluid2.append("CO2"); x0.append(0.99)  #cop=2.25
# =============================================================================

fluids_list=[fluid1+e for e in fluid2]
fluids=fluids_list[0]


n=200
eps=0.8
T_gc=40
#10**np.linspace(0,np.log10(100),10)    =   np.logspace(0,np.log10(100),10)
# =============================================================================
# P_gc=95#np.linspace(70,110,n)#95
# #P_gc=np.linspace(110,70,n)#95
# cop=np.zeros(n)
# x=np.linspace(0.85,0.9,n)
# =============================================================================
T_eva=-5
T_sep=10
T_ref=273.15
eta_c=0.8
xx=[]
copcop=[]


n=30
m=30
x=np.linspace(0.8,1,n)
P_gc=np.linspace(70,110,m)
y=P_gc#np.linspace(34,59,m)
cop=np.zeros([m,n])


glide=10

def process(P_gcc):
    funz=Funz(eps,P_gcc,T_gc,T_eva,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    print(T[8]-T_ref)
    return cop ,  T[8]-T[7]

def process7(P_gcc):
    funz=Funz7(eps,P_gcc,T_gc,T_eva-glide,T_sep,eta_c,mix,mix_l,mix_g)
        
    T,P,H,S,m,cop=funz.imp()
    T,P,H,S=funz.punti_rimanenti(T,P,H,S)
    print(T[8]-T_ref)
    return cop ,  T[8]-T[7]
    
def f_glide(x,args):

    P_gcc=args
    mix.set_mass_fractions([x, 1-x])
    cop,dT=process(P_gcc)
    print(x,dT)
    if dT<glide:
        cop,dT=process7(P_gcc)
        print(dT)
    return cop#,dT-glide
    
    

mix   = CP.AbstractState(lib, fluids)
mix_l = CP.AbstractState(lib, fluids)
mix_g = CP.AbstractState(lib, fluids)
    
# =============================================================================
# for i in range(n):
#     print('i =',i)
#     cop[i]=f_glide(x[i],P_gc)
#     
#   
# plt.figure(dpi=200)
# plt.plot(x,cop)
# plt.grid()
# =============================================================================
   
xx=[]
yy=[]
zz=[]
for i in range(n):
    print('\n\ni=',i,'\n'+'*'*80)
    for j in range(m):
        print('j=',j)
        try:
            if x[i]==1:
                mix   = CP.AbstractState(lib, 'CO2&CO2')
                mix_l = CP.AbstractState(lib, 'CO2&CO2')
                mix_g = CP.AbstractState(lib, 'CO2&CO2')
                copp=f_glide(0.99, y[j])
                cop[j,i] =copp
                xx.append(x[i])
                yy.append(y[j])
                zz.append(copp)
            else:
                copp=f_glide(x[i], y[j])
                cop[j,i] =copp
                xx.append(x[i])
                yy.append(y[j])
                zz.append(copp)
        except Exception as e:
            print(e)
z=cop
xi,yi = np.meshgrid(x,y)
metodo='linear'
zi = griddata((xx,yy),zz,(xi,yi),method=metodo)

# plot
plt.figure(dpi=200)
plt.contour(xi*100, yi, zi, levels=35, linewidths=0.5, colors='k')
cntr1 = plt.contourf(xi*100, yi, zi, levels=14, cmap="RdBu_r")
plt.colorbar(cntr1)
plt.xlabel(' x [%$kg_{CO_{2}}/kg_{tot}$]')#,fontsize=16)
plt.ylabel('$ P_{gc} $ [bar]')#,fontsize=16)
plt.title(' $ T_{gc} $=40°C, $ T_{eva} $='+str(T_eva-glide)+'°C, dT='+str(glide)+'°C,  '+fluids)  





plt.figure(dpi=200)
plt.contourf(xi*100,yi,zi, 60, cmap='jet')#,np.arange(0,1.01,0.01))
plt.colorbar();
#plt.plot(x,y,'k.')
plt.xlabel(' x [%$kg_{CO_{2}}/kg_{tot}$]')#,fontsize=16)
plt.ylabel('$ P_{gc} $ [bar]')#,fontsize=16)

plt.title(' $ T_{gc} $=40°C, $ T_{eva} $='+str(T_eva-glide)+'°C, dT='+str(glide)+'°C,  '+fluids)  


plt.figure(dpi=200)
plt.contourf(xi*100,yi,z, 30, cmap='jet')#,np.arange(0,1.01,0.01))
plt.colorbar();
#plt.plot(x,y,'k.')
plt.xlabel(' x [%$kg_{CO_{2}}/kg_{tot}$]')#,fontsize=16)
plt.ylabel('$ P_{gc} $ [bar]')#,fontsize=16)

plt.title('non interpolato')  

# =============================================================================
# fig, ax = plt.subplots()
# CS = ax.contour(xi, yi, zi)
# ax.clabel(CS, inline=1, fontsize=10)
# ax.set_title('Simplest default with labels')
# =============================================================================
