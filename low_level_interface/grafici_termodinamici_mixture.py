from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP

import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d

    

def curva_limitePH(mix, color, label, pc):
    n=20
    if pc==0:
        pc = mix.keyed_output(CP.iP_critical)
    #pc = 30*10**5#mix.keyed_output(CP.iP_critical)
    #tc = mix.keyed_output(CP.iT_critical)    

    pt=101325*2
    H_l = np.zeros(n)
    P_l = np.logspace(np.log10(pt), np.log10(pc),n)#np.linspace(pt,pc,n)  #np.logspace(np.log10(pt), np.log10(pc),n))
    for i in range(n):
        mix.update(CP.PQ_INPUTS, P_l[i], 0)
# =============================================================================
#         if i == n-1:
#             mix.update(CP.PT_INPUTS, pc, tc)
# =============================================================================
        H_l[i]=mix.hmass()
        
    plt.plot(H_l/1000,P_l/100000,color)    

    H_g = np.zeros(n)
    P_g = np.logspace(np.log10(pt), np.log10(pc),n)#np.linspace(pt,pc,n)  #list(np.logspace(np.log10(pt), np.log10(pc),n))
    for i in range(n):
        mix.update(CP.PQ_INPUTS, P_g[i], 1)
# =============================================================================
#         if i == n-1:
#             mix.update(CP.PT_INPUTS, pc, tc)
# =============================================================================
            

        H_g[i]=mix.hmass()

    #f = interp1d(P_g, H_g)#, kind='cubic')
    #P_g_n=np.linspace(pt,pc,100)
    plt.plot(H_g/1000,P_g/100000,color, label=label)   
    #plt.plot(f(P_g_n)/1000,P_g_n/100000,'k')   


def curva_limiteTS(mix, color, label, pc):
    n=20
    if pc==0:
        pc = mix.keyed_output(CP.iP_critical)
    #pc = 30*10**5#mix.keyed_output(CP.iP_critical)
    #tc = mix.keyed_output(CP.iT_critical)    

    pt=101325*2
    T_l = np.zeros(n)
    S_l = np.zeros(n)
    P_l = np.logspace(np.log10(pt), np.log10(pc),n)#np.linspace(pt,pc,n)  #np.logspace(np.log10(pt), np.log10(pc),n))
    for i in range(n):
        mix.update(CP.PQ_INPUTS, P_l[i], 0)
# =============================================================================
#         if i == n-1:
#             mix.update(CP.PT_INPUTS, pc, tc)
# =============================================================================
        #H_l[i]=mix.hmass()
        S_l[i]=mix.smass()
        T_l[i]=mix.T()
        
    plt.plot(S_l/1000,T_l,color)    

    T_g = np.zeros(n)
    S_g = np.zeros(n)
    P_g = np.logspace(np.log10(pt), np.log10(pc),n)#np.linspace(pt,pc,n)  #list(np.logspace(np.log10(pt), np.log10(pc),n))
    for i in range(n):
        mix.update(CP.PQ_INPUTS, P_g[i], 1)
# =============================================================================
#         if i == n-1:
#             mix.update(CP.PT_INPUTS, pc, tc)
# =============================================================================
            
        S_g[i]=mix.smass()
        T_g[i]=mix.T()
        #H_g[i]=mix.hmass()

    #f = interp1d(P_g, H_g)#, kind='cubic')
    #P_g_n=np.linspace(pt,pc,100)
    plt.plot(S_g/1000,T_g,color, label=label)   
    #plt.plot(f(P_g_n)/1000,P_g_n/100000,'k')   

def isoentropiche():
    len_s=11
    len_p=30
    S=np.linspace(1,2,len_s)*1000
    P=np.linspace(20,100,len_p)*100000
    H=np.zeros(len_p)
    for i in range (len_s):
        for j in range(len_p):
            H[j]=PropsSI('H','P',P[j],'S',S[i],'CO2')
        plt.plot(H/1000,P/100000,'b',linewidth=0.5)
        
def isoterme():
    passo=10
    len_p=30
    T_crit=PropsSI('Tcrit','CO2')
    T1=np.arange(T_crit,0+273.15,-passo)
    T2=np.arange(T_crit,180+273.15,passo)
    P=np.linspace(20,100,len_p)*100000
    H=np.zeros(len_p)
    for i in range (1,len(T2)):
        for j in range(len_p):
            H[j]=PropsSI('H','P',P[j],'T',T2[i],'CO2')
        plt.plot(H/1000,P/100000,'g',linewidth=0.5)
    
    len_p=200
    P=np.linspace(20,100,len_p)*100000
    H=np.zeros(len_p)
    for j in range(len_p):      
        H[j]=PropsSI('H','P',P[j],'T',T_crit,'CO2')
    plt.plot(H/1000,P/100000,'g',linewidth=0.5)
        
    for i in range (1,len(T1)):
        P_sat=PropsSI('P','Q',0,'T',T1[i],'CO2')
        len_p=15
        P1=np.linspace(20*10**5,P_sat-10,len_p)
        P2=np.linspace(P_sat+10,100*10**5,len_p)
        H=np.zeros(len_p)
        for j in range(len_p):
            H[j]=PropsSI('H','P',P1[j],'T',T1[i],'CO2')
        plt.plot(H/1000,P1/100000,'g',linewidth=0.5)
        
        H=np.zeros(len_p)
        for j in range(len_p):
            H[j]=PropsSI('H','P',P2[j],'T',T1[i],'CO2')
        plt.plot(H/1000,P2/100000,'g',linewidth=0.5)
    
        plt.plot(np.array([PropsSI('H','Q',0,'T',T1[i],'CO2'),PropsSI('H','Q',1,'T',T1[i],'CO2')])/1000,np.array([P_sat,P_sat])/100000,'g',linewidth=0.5)
       
        
        
def grafico_PH(P,H):
    plt.figure(dpi=200)
    curva_limitePH()
    isoentropiche()
    isoterme()
    
    plt.plot(H[:7],P[:7],'r')
    plt.plot([H[3],H[7],H[8],H[0]],[P[3],P[7],P[8],P[0]],'r')
    k=0
    xx=np.array([1,5,5,5,5,5,-14,-14,4])
    yy=np.array([1,-1,1,1,1,-4,1.5,-4,-4])
    for i,j in zip(H[:11],P[:11]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
        
        
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H')
    plt.grid()
    plt.show()  
    
def grafico_PH_sep(P,H,mix, mix_l, mix_g):
    plt.figure(dpi=200)
    pc = 30*10**5
    curva_limitePH(mix,'k','R455a',pc)
    curva_limitePH(mix_l, 'b','$ liquid_{fraction} $',pc)
    curva_limitePH(mix_g, 'orange','$ gas_{fraction} $',pc)
    #isoentropiche()
    #isoterme()
    
    plt.plot(H[:7],P[:7],'r')
    plt.plot([H[3],H[7],H[8],H[0]],[P[3],P[7],P[8],P[0]],'r')
    plt.plot([H[0],H[6]],[P[0],P[6]],'r')
    k=0
    xx=np.array([3,5,-10,5,5,-10,-10,3,3])
    yy=np.array([0.1,0,0,0.1,0.1,0.1,0.1,0,0])
    for i,j in zip(H[:11],P[:11]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
        
        
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H')
    plt.grid()
    plt.legend()
    plt.show()  
    
def grafico_PH_sep_IHX(P,H,mix, mix_l, mix_g,pc,fluids,x):
    plt.figure(dpi=200)
    #pc=65*10**5
    
    curva_limitePH(mix_l, 'b','$ liquid_{fraction} $',pc)
    curva_limitePH(mix_g, 'orange','$ gas_{fraction} $',pc)
    curva_limitePH(mix,'k','mix',pc)
    #isoentropiche()
    #isoterme()
    
# =============================================================================
#     plt.plot(H[:7],P[:7],'r')
#     plt.plot([H[3],H[7],H[8],H[0]],[P[3],P[7],P[8],P[0]],'r')
#     plt.plot([H[0],H[6]],[P[0],P[6]],'r')
#     k=0
#     xx=np.array([3,5,-10,5,5,-10,-10,3,3])
#     yy=np.array([0.1,0,0,0.1,0.1,0.1,0.1,0,0])
#     for i,j in zip(H[:11],P[:11]):
#         plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
#         plt.plot(i,j,'o',color='red')
#         k=k+1
# =============================================================================
        
    plt.plot(H[:9],P[:9],'r')
    plt.plot((H[5],H[9]),(P[5],P[9]),'r')
    plt.plot((H[9],H[10]),(P[9],P[10]),'r')
    k=0
    xx=np.array([1,5,5,5,5,5,-14,-14,4,5,-15])
    yy=np.array([1,-1,1,1,1,-4,1.5,-4,-4,1,-4])
    for i,j in zip(H[:11],P[:11]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H '+fluids+' mass fraction='+str(x))
    plt.grid()
    plt.legend()
    plt.show()  
    
def grafico_TS_sep_IHX(T,S,mix, mix_l, mix_g,pc,fluids,x):
    plt.figure(dpi=200)
    #pc=65*10**5
    P=T
    H=S
    
    curva_limiteTS(mix_l, 'b','$ liquid_{fraction} $',pc)
    curva_limiteTS(mix_g, 'orange','$ gas_{fraction} $',pc)
    curva_limiteTS(mix,'k','mix',pc)
 
        
    plt.plot(H[:9],P[:9],'r')
    plt.plot((H[5],H[9]),(P[5],P[9]),'r')
    plt.plot((H[0],H[10]),(P[0],P[10]),'r')
    plt.plot((H[9],H[10]),(P[9],P[10]),'r')
    k=0
    xx=np.array([0,0,0,0,0,0,0,0,0,0,0])
    yy=np.array([0,0,0,0,0,0,0,0,0,0,0])
    for i,j in zip(H[:11],P[:11]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
   
    plt.xlabel("S [kJ/kg]")
    plt.ylabel("T [K]")
    plt.title('Diagramma T-S '+fluids+' mass fraction='+str(x))
    plt.grid()
    plt.legend()
    plt.show()  
    
def grafico_PH_semplice(P,H):
    plt.figure(dpi=200)
    curva_limitePH()
    isoentropiche()
    isoterme()
    
    plt.plot(H[:6],P[:6],'r')
    plt.plot((H[0],H[5]),(P[0],P[5]),'r')
    
    k=0
    xx=np.array([1,5,5,5,5,5,-14,-14,4,5,3])
    yy=np.array([1,-1,1,1,1,-4,1.5,-4,-4,-4,1])
    for i,j in zip(H[:6],P[:6]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
        
        
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H')
    plt.grid()
    plt.show()  
 
    
def grafico_PH_sperimentale(P,H):
    plt.figure(dpi=200)
    curva_limitePH()
    isoentropiche()
    isoterme()
    
    plt.plot(H[:6],P[:6],'r')
    plt.plot((H[5],H[6]),(P[5],P[6]),'r')
    #plt.plot((H[4],H[6]),(P[4],P[6]),'r')
    plt.plot((H[4],H[6],H[7],H[8],H[9],H[0]),(P[4],P[6],P[7],P[8],P[9],P[0]),'r')
    #plt.plot(H[6:10],P[6:10],'r')
    
    k=0
    xx=np.array([1,5,5,5,5,-10,0,2,-10,-10])
    yy=np.array([1,-1,1,1,1,-5,-5,-5,-5,1])
    for i,j in zip(H[:10],P[:10]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
        
        
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H')
    plt.grid()
    plt.show()  
    
    
def grafico_PH_eiettore_sep(P,H,mix, mix_l, mix_g,pc,fluids,x):
    plt.figure(dpi=200)
    #pc=65*10**5
    
    curva_limitePH(mix_l, 'b','$ liquid_{fraction} $',pc)
    curva_limitePH(mix_g, 'orange','$ gas_{fraction} $',pc)
    curva_limitePH(mix,'k','mix',pc)
    #isoentropiche()
    #isoterme()

        
    plt.plot(H[:6],P[:6],'r')
    plt.plot((H[0],H[6],H[10]),(P[0],P[6],P[10]),'r')
    plt.plot(H[6:10],P[6:10],'r')
    plt.plot((H[9],H[10]),(P[9],P[10]),'r')
    k=0
    xx=np.array([1,5,5,5,5,5,-14,-14,4,5,3])
    yy=np.array([1,-1,1,1,1,-4,1.5,-4,-4,-4,1])
    for i,j in zip(H[:11],P[:11]):
        plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color='red')
        k=k+1
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H '+fluids+' mass fraction='+str(x))
    plt.grid()
    plt.legend()
    plt.show()  