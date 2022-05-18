#poison.py
import numpy as np
import matplotlib.pyplot as plt

#Constants
gamma = {'I':0.0639,'Xe':0.00237}
lamda = {'I':2.87e-5, 'Xe':2.09e-5} #sec^1
phi = 1.1e14 #n/cm2/s
#phi = 1.1e14 #n/cm2/s

sigma_F = 587e-24 #cm2
N = 0.04833e24 #atoms/cm3
Sigma_F = sigma_F*N

sigma_Xe = 2.6e-18


#Define Functions
def dI_dt(I):
    return gamma['I']*Sigma_F*phi-lamda['I']*I

def dXe_dt(I,Xe):
    return lamda['I']*I+gamma['Xe']*Sigma_F*phi-lamda['Xe']*Xe-Xe*sigma_Xe*phi

dt=1 #sec

#Loop
begin = 0 #hours
end = 500 #hours
time = np.arange(begin,end*3600+1) #sec

for t in time:
    if t == 200*3600:
        phi = 0
        
    if t==0: #Iniitial Conditions for ODEs
        I = [0]
        Xe= [0]
    else:
        i=I[-1]+dI_dt(I[-1])*dt
        xe=Xe[-1]+dXe_dt(I[-1],Xe[-1])*dt
        I.append(i)
        Xe.append(xe)

ax =  plt.gca()
plt.xlabel('time since step change, t (hr)')
plt.ylabel('Concentration (atoms/cubic centimeter)')
plt.yscale('log')
plt.ylim(1e16,1e20)
plt.xlim(-1,500)
#plt.yticks(np.array([1e18,5e18,1e19]))
plt.plot(time/3600,Xe,label='Xenon 135')
plt.plot(time/3600,I,label='Iodine 135')
plt.grid()

plt.legend(loc='best')
twin_ax = ax.twinx()
twin_ax.set_yscale('log')
twin_ax.set_ylim(ax.get_ylim())



plt.show()
