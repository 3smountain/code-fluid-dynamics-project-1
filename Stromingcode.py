# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 16:41:44 2023

@author: dries
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import time
import cmath

#Voor opdracht a.
#Dit is de funcite die circels maakt.
def circlegen(N,z0,R):
    arr=np.linspace(0,2*np.pi,N)
    z1=np.exp(arr*1j)
    z2=R*z1+z0
    return z2
#Dit is de functie die de joukoutransformatie uitvoert.
def joukou(z):
    s=z+1/z
    return s
#Hier maken we de circel met het middelpunt -0.1+0.22i en de straal 1.12.
circle=circlegen(1000,-0.1+0.22j,1.12)
plt.figure()
plt.plot(np.real(circle),np.imag(circle), markersize=0.5)
plt.axis('square')
plt.show()
#Hier maken we de wingshape door de circel te transformeren.
wing=joukou(circle)
plt.figure()
plt.plot(np.real(wing),np.imag(wing), 'g.',markersize=0.5)
plt.axis([-3,3, -2, 2])
plt.title('joukowski transformation of cylinder')
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.savefig('joukoplanewing.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
#%%
#Voor opdracht b.
#Hier definieren we een functie die een grid maakt rondom de circel.
def gridgen(Nj, Nk, R):
    R=np.linspace(1.12,R,Nj)
    ang=np.exp(1j*np.linspace(0,2*np.pi,Nk))
    xv, yv = np.meshgrid(R,ang)
    arr = xv*yv-0.1+0.22j
    return arr 
#Hier definieren we een functie die het potentiaal berekend.
def comppot(z,gamma):
    O=0.1-0.22j
    fz= z+O+(1.12**2)/(z+O)-gamma*1j*np.log(z+O)/(2*np.pi)
    return fz 
#Hier maken we het grid aan 
GRID=gridgen(1500,1500,3)
a = np.real(GRID)
b = np.imag(GRID)
c= np.imag(comppot(GRID,-3))
d= np.imag(comppot(GRID,0))
levels = np.linspace(-4,4,15)
plt.figure()
#Hier plotten we de stroomlijnen voor gamma=-3 en gamma=0
plt.contour(a,b,c, levels)
plt.title('Imaginary part of velocity potential for cylinder with Gamma=-3')
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(circle),np.imag(circle),'r.', markersize=0.5)
plt.savefig('cylinderwithgamma-3.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.figure()
plt.title('Imaginary part of velocity potential for cylinder with Gamma=0')
plt.contour(a,b,d, levels)
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(circle),np.imag(circle),'r.', markersize=0.5)
plt.savefig('cylinderwithgamma0.pdf', dpi=4000, format='pdf',bbox_inches='tight')

#%%
#Voor opdracht c.
#Hier plotten we de stroomlijnen, maar dan voor de vleugel.
e0 = joukou(GRID)
f0 = joukou(GRID)
e=np.real(e0)
f=np.imag(f0)
plt.figure()
plt.plot(e,f)
plt.show()
g= np.imag(comppot(GRID,-3))
h= np.imag(comppot(GRID,0))
levels = 15
plt.figure()
plt.contour(e,f,g, levels)
plt.title('Imaginary part of velocity potential for wingshape with Gamma=-3')
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(wing),np.imag(wing),'r.', markersize=0.5)
plt.savefig('wingwithgamma-3.pdf', dpi=4000, format='pdf',bbox_inches='tight')

plt.figure()
plt.contour(e,f,h, levels)
plt.title('Imaginary part of velocity potential for wingshape with Gamma=0')
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(wing),np.imag(wing),'r.', markersize=0.5)
plt.savefig('wingwithgamma0.pdf', dpi=4000, format='pdf',bbox_inches='tight')
#%%
#Voor opdracht d.
#Hier plotten we de stroomlijnen bij de scherpe punt van de vleugel.
levels=25
test3=np.imag(comppot(GRID,3))
test0=np.imag(comppot(GRID,0))
testmin3=np.imag(comppot(GRID,-3))
plt.figure()
plt.contour(e,f,testmin3, levels)
plt.title('Imaginary part of velocity potential at sharp edge with Gamma=-3')
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(wing),np.imag(wing),'r.', markersize=0.5)
plt.axis([1.5,2.5,-1,1])
plt.savefig('edgewithgamma-3.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.figure()
plt.contour(e,f,test0, levels)
plt.title('Imaginary part of velocity potential at sharp edge with Gamma=0')
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(wing),np.imag(wing),'r.', markersize=0.5)
plt.axis([1.5,2.5,-1,1])
plt.savefig('edgewithgamma0.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.figure()
plt.contour(e,f,test3, levels)
plt.title('Imaginary part of velocity potential at sharp edge with Gamma=3')
plt.colorbar()
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(wing),np.imag(wing),'r.', markersize=0.5)
plt.axis([1.5,2.5,-1,1])
plt.savefig('edgewithgamma3.pdf', dpi=4000, format='pdf',bbox_inches='tight')
#%%
#Voor opdracht e.
#Hier plotten we weer de stroomlijnen voor gamma=-3, maar deze keer ook het reeele deel van de potentiaal.
levelse=100
teste0=comppot(GRID,-3)
teste1=np.imag(teste0)
teste2=np.real(teste0)
plt.figure()
plt.contour(e,f,teste1, levelse)
plt.colorbar()
plt.contour(e,f,teste2,levelse)
plt.colorbar()
plt.title('both real and imaginary velocity potential at sharp edge with Gamma=-3')
plt.xlabel('Real part of position')
plt.ylabel('Imaginary part of position')
plt.plot(np.real(wing),np.imag(wing),'r.', markersize=0.5)
plt.axis([1.5,2.5,-1,1])
plt.savefig('realandimagedgewithgamma-3.pdf', dpi=4000, format='pdf',bbox_inches='tight')
#%%
#voor f.
#Hier bepalen we het snelheidsveld.
def vel(grid,pot):
    n=len(grid)
    w=np.zeros((n,n),dtype=np.complex128)
    for i in range(n):
        w[i,-1]=(pot[i,-1]-pot[i,-2])/(grid[i,-1]-grid[i,-2])
    for i in range(n):
        for j in range(n-1):
            w[i,j]=(pot[i,j+1]-pot[i,j])/(grid[i,j+1]-grid[i,j])
    v_x=w.real
    v_y=-w.imag
    return v_x,v_y
levels=25
poten=comppot(GRID,-3)
gridjou=joukou(GRID)
vx,vy=vel(GRID,poten)
#Hier plotten we het snelheidsveld nog even.
plt.figure()
plt.contourf(np.real(GRID),np.imag(GRID),vx, levels)
plt.colorbar()
plt.figure()
plt.contourf(np.real(GRID),np.imag(GRID),vy, levels)
plt.colorbar()
#Hier bepalen en plotten we de druk van de circel.
p=-(vx**2+vy**2)/2
plt.figure()
plt.contourf(np.real(GRID),np.imag(GRID),p, levels)
plt.colorbar()
plt.savefig('pressureofcylwithgamma-3.pdf', dpi=4000, format='pdf',bbox_inches='tight')
#Hier moeten we de levels wat aanpassen om uitschieters niet te plotten.
levels0=np.arange(-1,2,0.01)
levels1=np.arange(-1.5,0,0.01)
#Hier plotten we het snelheidsveld van de wing.
vx1,vy1=vel(gridjou,poten)
plt.figure()
plt.contourf(np.real(gridjou),np.imag(gridjou),vx1, levels0)
plt.colorbar()
plt.figure()
plt.contourf(np.real(gridjou),np.imag(gridjou),vy1, levels0)
plt.colorbar()
#Hier bepalen en plotten we de druk van de wing.
p1=-(vx1**2+vy1**2)/2
plt.figure()
plt.contourf(np.real(gridjou),np.imag(gridjou),p1, levels1)
plt.colorbar()
plt.savefig('pressureofwingwithgamma-3.pdf', dpi=4000, format='pdf',bbox_inches='tight')

#%% 
#Voor g.
#Hier integreren we over de rand van een vorm (eerst circel dan de vleugel) en bepalen we daarmee de kracht.
#Eigenlijk someren we in plaats van integreren.
circlepot=comppot(circle,-3)
def linint(line,pot):
    n=len(line)
    w=np.zeros(n,dtype=np.complex128)
    w[-1]=(pot[-1]-pot[-2])/(line[-1]-line[-2])
    for i in range(n-1):
        w[i]=(pot[i+1]-pot[i])/(line[i+1]-line[i])
    somel=np.zeros(n,dtype=np.complex128)
    somel[0]=w[0]**2*(line[1]-line[-1])/2
    for i in range(n-2):
        somel[i+1]=w[i+1]**2*(line[i+2]-line[i])/2
    somel[-1]=w[-1]**2*(line[0]-line[-2])/2
    somelabs=np.absolute(somel)
    print(np.mean(somel))
    #Hier verwijderen we wat punten die uitschieters zijn.
    for i in range(n):
        if np.absolute(somel[i])>np.mean(somelabs)+np.std(somelabs):
            somel[i]=0
    som=1j/2*np.sum(somel)
    return som
Fcircle=linint(circle,circlepot)
Fwing=linint(wing,circlepot)
print(Fcircle,Fwing)