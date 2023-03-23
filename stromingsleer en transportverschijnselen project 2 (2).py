#!/usr/bin/env python
# coding: utf-8

# In[93]:


#Author: Dries van de Bult
import matplotlib.animation
import numpy as np
import matplotlib.pyplot as plt


# In[94]:


#Here we define the constants for the analytical part of the project.
U0=1
T=10
rho=1
mu=1
omega=2*np.pi/T
k=np.sqrt(omega*rho/(2*mu))
def tau(z):
    temp=(1/omega)*k*z
    return temp


# In[95]:


#code inspired by "Kolibril"'s answer on https://stackoverflow.com/questions/43445103/inline-animations-in-jupyter
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams['figure.dpi'] = 150  
plt.ioff()
fig, ax = plt.subplots()

#Here we make the animated plot for exercise e of the analytical part of project 2
z= np.linspace(0,5,100)

def animate(t):
    plt.cla()
    u=np.exp(-k*z)*np.sin(2*np.pi*(t-tau(z))/T)
    plt.plot(z,u)
    plt.xlim(0,5)
    plt.ylim(-1,1)

matplotlib.animation.FuncAnimation(fig, animate, frames=100)


# In[96]:


#Here we also make a contour plot for the analytical part of the project.
def u(z,t):
    u=np.zeros((len(z),len(t)))
    for i in range(len(t)):
        u[i]=np.exp(-k*z)*np.sin(2*np.pi*(t[i]-tau(z))/T)
    return u
z=np.linspace(0,5,100)
t=np.linspace(0,10,100)
levels=25
uzt0=u(z,t)
plt.figure()
plt.contourf(z,t,uzt0,levels)
plt.xlabel('z')
plt.ylabel('t')
plt.title('controur plot for e of the analytical part')
plt.savefig('controurofeanalyticalpart.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[97]:


#Here we define the needed constants for the numerical part of the project.
rho=1000
mu=1
T=0.1
omega=2*np.pi/T
U0=1
dz=0.001#This is our choise for dz
H=0.03#This is our choise for H
Nz=int(H/dz)

k=np.sqrt(rho*omega/(2*rho))
dt=1/10*1/(mu/rho)*dz**2#by defining dt like this we guarantee that dt is smaller then dz
Nt=10*T/dt#We've chosen to plot ten times the period of the function
Nt=int(Nt)


# In[98]:


#Here we define the second derivative approximation function
def dzz(u):
    dzzu=(np.roll(u,-1)+np.roll(u,1)-2*u)/dz**2
    return dzzu
#Here we define the function that calculates u over time. In this case the lower plate moves as U0sin(omega t)
def timestep(uu):
    for t in range(Nt):
        u=uu[t]
        F=(mu/rho)*dzz(u)
        uu[t+1]=uu[t]+dt*F
        uu[t+1][Nz]=0
        uu[t+1][0]=U0*np.sin(omega*t*dt)
    return uu


# In[99]:


#Here we run the functions of above for our chosen Nt and Nz.
uzt=np.zeros((Nt+1,Nz+1))
uzt1=timestep(uzt)

#Here we define the axis we will use
time=np.arange(0,Nt+1)*dt
zaxis=np.arange(0,Nz+1,1)*dz
levels=25


# In[100]:


#Here we plot the contour of u agianst z and t
#This plot is a bit weird sometimes. I've had multiple times that it showed another plot as well as the contour plot or sometimes another plot instead of the right one.
#Because of the above comment you may want to run this cell a couple of times. With me it works correctly after the thirth time running.
plt.figure()
plt.contourf(zaxis,time,uzt1,levels)
plt.ylabel('time')
plt.xlabel('z')
plt.colorbar()
plt.title('contour plot of u against z and time with the 1st plate motion')
plt.savefig('controuroffirstmotion.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[101]:



plt.figure()
plt.plot(zaxis,uzt1[0],'r.')
plt.ylabel('u')
plt.xlabel('z')
plt.title('u against z for t=0 (before movement lower plate) with the 1st plate motion')
plt.savefig('1stplateuvsztis0.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(zaxis,uzt1[int(Nt/30)],'r.')
plt.ylabel('u')
plt.xlabel('z')
plt.title('u against z for t=Nt/30 (at positive peak) with the 1st plate motion')
plt.savefig('1stplateuvsztisNtdiv30.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(zaxis,uzt1[-1],'r.')
plt.title('u against z for t=Nt (in between negative and positive peak) with the 1st plate motion')
plt.ylabel('u')
plt.xlabel('z')
plt.savefig('1stplateuvsztisNt.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[102]:


plt.figure()
plt.plot(time,uzt1[:,0], 'b.')
plt.ylabel('u')
plt.xlabel('time')
plt.title('u against time for z=0 with the 1st plate motion')
plt.savefig('1stplateuvstzis0.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(time,uzt1[:,int(Nz/2)], 'b.')
plt.ylabel('u')
plt.xlabel('time')
plt.title('u agianst time for z=Nz/2 (halfway between the two plates) with the 1st plate motion')
plt.savefig('1stplateuvstzisNzdiv2.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(time,uzt1[:,(Nz-1)], 'b.')
plt.ylabel('u')
plt.xlabel('time')
plt.title('u against time for z=Nz-1 (almost at the upper plate) with the 1st plate motion')
plt.savefig('1stplateuvstzisNzmin1.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[103]:



#%%
#Here we define the timestep function for a plate that alternates between moving at U0 and -U0 velocity.
T2=T/2
def timestep2(uu):
    for t in range(Nt):
        u=uu[t]
        F=(mu/rho)*dzz(u)
        uu[t+1]=uu[t]+dt*F
        uu[t+1][Nz]=0
        m=t*dt%T
        if m<T2:
            uu[t+1][0]=U0
        if m>T2:
            uu[t+1][0]=-U0
    return uu


# In[104]:


#Here we run the our second timestep function.
uzt=np.zeros((Nt+1,Nz+1))
uzt2=timestep2(uzt)
#here we define the axis we will use.
time=np.arange(0,Nt+1)*dt
zaxis=np.arange(0,Nz+1,1)*dz
levels=25


# In[105]:


#Here we plot the contour of u against z and t for the second timestep function.
plt.figure()
plt.contourf(zaxis,time,uzt2,levels)
plt.ylabel('time')
plt.xlabel('z')
plt.colorbar()
plt.title('contour plot of u against z and time with the 2nd plate motion')
plt.savefig('controurofsecondmotion.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[106]:


plt.figure()
plt.plot(zaxis,uzt2[0],'r.')
plt.ylabel('u')
plt.xlabel('z')
plt.title('u against z for t=0 (before movement lower plate) with the 2nd plate motion')
plt.savefig('2ndplateuvsztis0.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(zaxis,uzt2[int(Nt/30)],'r.')
plt.ylabel('u')
plt.xlabel('z')
plt.title('u against z for t=Nt/30 (at positive peak) with the 2nd plate motion')
plt.savefig('2ndplateuvsztisNtdiv30.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(zaxis,uzt2[-1],'r.')
plt.title('u against z for t=Nt (in between negative and positive peak) with the 2nd plate motion')
plt.ylabel('u')
plt.xlabel('z')
plt.savefig('2ndplateuvsztisNt.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[107]:


plt.figure()
plt.plot(time,uzt2[:,0], 'b.')
plt.ylabel('u')
plt.xlabel('time')
plt.title('u against time for z=0 with the 2nd plate motion')
plt.savefig('2ndplateuvstzis0.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(time,uzt2[:,int(Nz/2)], 'b.')
plt.ylabel('u')
plt.xlabel('time')
plt.title('u agianst time for z=Nz/2 (halfway between the two plates) with the 2nd plate motion')
plt.savefig('2ndplateuvstzisNzdiv2.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()
plt.figure()
plt.plot(time,uzt2[:,(Nz-1)], 'b.')
plt.ylabel('u')
plt.xlabel('time')
plt.title('u against time for z=Nz-1 (almost at the upper plate) with the 2nd plate motion')
plt.savefig('2ndplateuvstzisNzmin1.pdf', dpi=4000, format='pdf',bbox_inches='tight')
plt.show()


# In[ ]:




