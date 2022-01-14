from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

def curva_limitePH():
    T_lim=np.linspace(250,PropsSI('Tcrit','CO2'),400)
    l=len(T_lim)
    T_l=np.zeros(2*l)
    H_l=np.zeros(2*l)
    
    for i in range(l):
        H_l[i]=PropsSI('H','Q',0,'T',T_lim[i],'CO2')
        T_l[i]=T_lim[i]
    for j in range(l):
        H_l[j+l]=PropsSI('H','Q',1,'T',T_lim[l-j-1],'CO2')
        T_l[j+l]=T_lim[l-j-1]
            
    P_l=np.zeros(2*l)
        
    for i in range(l):
        P_l[i]=PropsSI('P','Q',0,'T',T_lim[i],'CO2')
        T_l[i]=T_lim[i]
    for j in range(l):
        P_l[j+l]=PropsSI('P','Q',1,'T',T_lim[l-j-1],'CO2')
        T_l[j+l]=T_lim[l-j-1]
    
    plt.plot(H_l/1000,P_l/100000,'k')    

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
    plt.title('Diagramma P-H')
    plt.grid()
    plt.show()  
    
def grafico_PH_sep(P,H):
    plt.figure(dpi=200)
    curva_limitePH()
    isoentropiche()
    isoterme()
    
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
    plt.title('Diagramma P-H')
    plt.grid()
    plt.show()  
    
def grafico_PH_semplice(P,H,color,number):
    curva_limitePH()
    isoentropiche()
    isoterme()
    
    plt.plot(H,P,'--',color=color)
    plt.plot((H[0],H[5]),(P[0],P[5]),'--',color=color)
    
    k=0
    xx=np.array([1,5,5,5,5,5,-14,-14,4,5,3])
    yy=np.array([1,-1,1,1,1,-4,1.5,-4,-4,-4,1])
    for i,j in zip(H,P):
        if number==1:
            plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color=color)
        k=k+1
        
        
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    plt.title('Diagramma P-H')
    plt.grid()
   
    
def grafico_PH_super_semplice(P,H,color,number):
    curva_limitePH()
    isoentropiche()
    isoterme()
    
    plt.plot(H,P,'--',color=color)
    plt.plot((H[0],H[3]),(P[0],P[3]),'--',color=color)
    
    
    k=0
    xx=np.array([5,5,5,5,5,5,-14,-14,4,5,3])
    yy=np.array([1,-1,1,1,1,-4,1.5,-4,-4,-4,1])
    for i,j in zip(H,P):
        if number==1:
            plt.annotate(str(k),xy=(i+xx[k],j+yy[k]))
        plt.plot(i,j,'o',color=color)
        k=k+1
        
   
    plt.xlabel("H [kJ/kg]")
    plt.ylabel("P [bar]")
    #plt.title('Diagramma P-H')
    plt.grid()
    
    
 
    
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