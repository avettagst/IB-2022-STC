from cProfile import label
from cmath import sqrt
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import lib
import math as m

gamma = 1.2 # (watt por km)^(-1)
P0 = 0.1 #potencia pico en W
alfa = 0.00414

Lef = 1/alfa


n = 2048
dt = 1 #paso de tiempo
t = np.linspace(-n*dt/2,n*dt/2,n-1)

f_s = 1/dt
f = np.linspace(-f_s/2,f_s/2,n-1)
f = sc.fft.fftshift(f)
f = np.roll(f,1)

t0 = 12.5
A0t = np.sqrt(P0)*np.exp(-0.5*(t/t0)**2)
A0f = sc.fft.fft(A0t)

L = 200 #largo de la fibra en km
d = 10 #cant de pasos para el vector z
z = np.linspace(0, L, d)

#no tenemos atenuación, por lo que Leff = L
for pos in z:
    Azt = A0t*np.exp(1j*gamma*pos*np.abs(A0t*np.exp(-0.5*alfa*pos))**2)*np.exp(-alfa*pos)
    Azf = sc.fft.fft(Azt)

phimax = gamma*P0*Lef
M = phimax/m.pi+0.5

print(phimax)
print(M)


plt.plot(t,np.abs(Azt)**2, 'r', label = 'Pulso en z = 500km')
plt.plot(t,np.abs(A0t)**2, 'k', label = 'Pulso en z = 0')
plt.title('Dominio del tiempo')
plt.xlabel('t[ps]')
plt.ylabel('|A(z, t)|^2')
plt.legend()
plt.show()

plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(A0f)**2), 'k--', lw = 0.5, label = 'Espectro inicial')
plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(Azf)**2), 'r', label = 'Espectro en z = L')
plt.title('Dominio de la frecuencia')
plt.xlabel('f[GHz]')
plt.ylabel('|A(z, f)|^2')
plt.xlim(-300,300)
#plt.ylim(0,0.1)
plt.title('L = ' + str(L) + ' km // P0 = ' + str(1000*P0) + ' mW')
plt.legend()
plt.show()

