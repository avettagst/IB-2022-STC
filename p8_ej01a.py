from cProfile import label
from cmath import sqrt
from turtle import color
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import lib
import math as m

#SPECS
alfa = 0.0414 #kilometro^(-1)
#b2 = -20.2 #picosegundos cuadrado sobre kilómetro
#b3 = 0.163 #ps^3/km
#gamma = 1.85 # (watt por km)^(-1)

#alfa = 0
b2 = 0
b3 = 0
gamma = 0




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
L = 50
p = 500
z = np.linspace(0,L,p+1)
d = z[1]-z[0]

print('Paso: ' +str(d)+' km') #debería ser L/p

#PULSO INICIAL
t0 = 10 #en picosegundos
A0t = np.exp(-0.5*(t/t0)**2) #pulso inicial
A0f = sc.fft.fft(A0t) #su espectro



#copio el pulso inicial porque en el for voy a ir reescribiendo el paso anterior (en un solo dominio)
Azf = A0f

#PARA VERIFICAR ATENUACIÓN
E0 = np.trapz(abs(A0f)**2,dx=f[1]-f[0])
energia = []
energia_teo = []


for pos in z:
    Azf = Azf*np.exp((1/6)*1j*b3*d*(2*np.pi*f)**3)*np.exp(0.5*1j*b2*d*(2*np.pi*f)**2)
    Azt = sc.fft.ifft(Azf)
    Azt = Azt*np.exp(1j*gamma*d*np.abs(Azt)**2)*np.exp(-0.5*d*alfa)
    Azf = sc.fft.fft(Azt)
    energia = np.append(energia, np.trapz(abs(Azf)**2,dx=f[1]-f[0]))
    energia_teo = np.append(energia_teo, E0*np.exp(-alfa*pos))






plt.plot(t,np.abs(Azt)**2, 'r')
plt.plot(t,np.abs(A0t)**2, 'k')
plt.title('Dominio del tiempo')
plt.xlabel('t[ps]')
plt.ylabel('|A(z=0, t)|^2')
#plt.legend()
plt.show()

plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(A0f)**2), 'k') 
plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(Azf)**2), 'r')
plt.title('Dominio de la frecuencia')
plt.xlabel('f[GHz]')
plt.ylabel('|A(z=0, w)|^2')
plt.xlim(-100,100)
plt.show()


plt.plot(z, 10*np.log10(energia/E0), 'bo')
plt.plot(z, 10*np.log10(energia_teo/E0), 'r')
plt.xlabel('z [km]')
plt.ylabel('Decaimiento de energía [dB]')
plt.show()

#Pendiente de la recta en dB (debe ser -alfa_dB)
print('Pendiente en Db:')
print(10*(np.log10(energia[-1])- np.log10(energia[0]))/(z[-1]))