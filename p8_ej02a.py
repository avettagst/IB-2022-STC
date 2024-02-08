from cProfile import label
from cmath import sqrt
from turtle import color
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import lib
import math as m

#SPECS
alfa = 0
#b2 = 0
#b3 = 0
gamma = 0

#alfa = 0.0414 #kilometro^(-1)
b2 = -20.2 #picosegundos cuadrado sobre kilómetro
b3 = 0.163 #ps^3/km
#gamma = 1.85 # (watt por km)^(-1)


#CANT DE MUESTRAS
n = 2048

#VECTOR DE TIEMPOS (ps)
dt = 0.5
t = np.linspace(-n*dt/2,n*dt/2,n-1)

#VECTOR DE FRECUENCIAS (THz)
fs = 1/dt
df = fs/n
f = np.linspace(-fs/2,fs/2,n-1)
f = sc.fft.fftshift(f)
f = np.roll(f,1)

#VECTOR DE DISTANCIAS (km)
L = 80
p = 500
z = np.linspace(0,L,p+1)
d = z[1]-z[0]

#PULSO INICIAL
t0 = 25 #en picosegundos
m = 2 #orden
A0t = np.exp(-0.5*(t/t0)**(2*m)) #pulso inicial
A0f = sc.fft.fft(A0t) #su espectro

#copio el pulso inicial porque en el for voy a ir reescribiendo el paso anterior (en un solo dominio)
Azf = A0f

for pos in z:
    Azf = Azf*np.exp((1/6)*1j*b3*d*(2*np.pi*f)**3)*np.exp(0.5*1j*b2*d*(2*np.pi*f)**2)
    Azt = sc.fft.ifft(Azf)
    Azt = Azt*np.exp(1j*gamma*d*np.abs(Azt)**2)*np.exp(-0.5*d*alfa)
    Azf = sc.fft.fft(Azt)


plt.plot(t,np.abs(Azt)**2, 'r', label = 'Pulso en z = L')
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
plt.title('L = ' + str(L) + ' km')
plt.legend()
plt.show()