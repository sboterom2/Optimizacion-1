from sympy import symbols 
import sympy as smp
from sympy import *
t, g = smp.symbols('t g')
m = smp.symbols('m:2')
L = smp.symbols('L:2')
z=smp.symbols('z:2')
qq=smp.symbols('q:4', cls=smp.Function)
the = smp.symbols('theta:2', cls=smp.Function)
n=2
thef=zeros(1,n)
thef_d=zeros(1,n)
thef_dd=zeros(1,n)
xn=zeros(1,n)
yn=zeros(1,n)
k=zeros(1,n)
p=zeros(1,n)
eel=zeros(1,n)
qqf=zeros(1,2*n)
qqfd=zeros(1,2*n)
for i in range(0, n):  
    #funciones
    temp=the[i]
    temp=temp(t)
    thef[0,i]=temp
    #derivada
    temp=thef[0,i]
    temp=smp.diff(temp, t)
    thef_d[0,i]=temp
    #segunda deriva
    temp=thef_d[0,i]
    temp=smp.diff(temp, t)
    thef_dd[0,i]=temp
    #x de la posicion
    temp=L[i] *smp.sin(thef[0,i]) 
    xn[0,i]=temp
    xn[0,i]=sum(xn)
    # y de la posicion 
    temp=L[i] *smp.cos(thef[0,i]) 
    yn[0,i]=temp
    yn[0,i]=sum(yn)       
yn=yn*-1 
for i in range(0, 2*n):
    #funciones estado
    temp=qq[i]
    temp=temp(t)
    qqf[0,i]=temp
    #derivada funciones de estado
    temp=qqf[0,i]
    temp=smp.diff(temp, t)
    qqfd[0,i]=temp 
    
for i in range(0, n):
    #cinetica
    temp=1/2*m[i]*(smp.diff(xn[0,i], t)**2 + smp.diff(yn[0,i], t)**2)
    k[0,i]=temp
    #potencial
    temp=m[i]*g*yn[0,i]
    p[0,i]=temp
kt=sum(k)
pt=sum(p)
la=kt-pt
for i in range(0,n):
    #ecuaciones euler-lagrange
    temp=smp.diff(smp.diff(la, thef_d[0,i]), t)-smp.diff(la, thef[0,i]).simplify()
    eel[0,i]=temp
#primera ecuacion    
a=eel[0,0].simplify()
   
f=a.replace(smp.sin(thef[0,0]), thef[0,0])
a=f.replace(smp.sin(thef[0,1]), thef[0,1])
f=a.replace(smp.cos(thef[0,0]), 1)
a=f.replace(smp.sin(thef[0,0]-thef[0,1]), 1)
f=a.replace(smp.cos(thef[0,1]), 1)
a=f.replace(smp.cos(thef[0,0]-thef[0,1]), 1)
a=a.simplify()

f=a.replace(thef_dd[0,0],qqfd[0,1])
a=f.replace(thef_d[0,0],qqf[0,1])
f=a.replace(thef[0,0],qqf[0,0])

a=f.replace(thef_dd[0,1],qqfd[0,3])
f=a.replace(thef_d[0,1],qqf[0,3])
a=f.replace(thef[0,1],qqf[0,2])

sol = smp.solve(a, (qqfd[0,1]), simplify=False, rational=False)
print(sol)
#segunda ecuacion











     
    
   

    
    







