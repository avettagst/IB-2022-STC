from cProfile import label
from cmath import sqrt
from turtle import color
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import lib
import math as m

#SPECS
#alfa = 0.0414 #kilometro^(-1)
b2 = 10 #picosegundos cuadrado sobre kilómetro
#b3 = 0.163 #ps^3/km
gamma = 1.8 #(watt por km)^(-1)

alfa = 0
#b2 = 0
b3 = 0
#gamma = 0




#CANT DE MUESTRAS
n = 8192

#VECTOR DE TIEMPOS (ps)
dt = 0.5
t = np.linspace(-n*dt/2,n*dt/2,n-1)

#VECTOR DE FRECUENCIAS (THz)
fs = 1/dt
#df = fs/n
f = np.linspace(-fs/2,fs/2,n-1)
f = sc.fft.fftshift(f)
f = np.roll(f,1)

#VECTOR DE DISTANCIAS (km)
L = 100
p = 1000
z = np.linspace(0,L,p+1)
d = z[1]-z[0]

print('Paso: ' +str(d)+' km') #debería ser L/p

#PULSO INICIAL
t0 = 30 #picosegundos
P0 = 0.05 #watts. es la potencia de cada pulso 
df = 0.05 #THz
A01 = sqrt(P0)*(np.exp(-0.5*(t/t0)**2))
A02 = sqrt(P0)*(np.exp(-0.5*(t/t0)**2))*np.exp(1j*2*m.pi*df*t)
A0t = A01 + A02
A01f = A0f = sc.fft.fft(A01)
A0f = sc.fft.fft(A0t) #su espectro



#copio el pulso inicial porque en el for voy a ir reescribiendo el paso anterior (en un solo dominio)
Azf = A0f
Azf1 = A01f

for pos in z:
    Azf = Azf*np.exp((1/6)*1j*b3*d*(2*np.pi*f)**3)*np.exp(0.5*1j*b2*d*(2*np.pi*f)**2)
    Azt = sc.fft.ifft(Azf)
    Azt = Azt*np.exp(1j*gamma*d*np.abs(Azt)**2)*np.exp(-0.5*d*alfa)
    Azf = sc.fft.fft(Azt)

    Azf1 = Azf1*np.exp((1/6)*1j*b3*d*(2*np.pi*f)**3)*np.exp(0.5*1j*b2*d*(2*np.pi*f)**2)
    Azt1 = sc.fft.ifft(Azf1)
    Azt1 = Azt1*np.exp(1j*gamma*d*np.abs(Azt1)**2)*np.exp(-0.5*d*alfa)
    Azf1 = sc.fft.fft(Azt1)

#Calculo la diferencia de altura entre el pico del canal y el pico más cercano
e1 = 1
for i in range(len(f)):
    if(abs(f[i] - df) < e1):
        e1 = abs(f[i] - df)
        pos1 = i

e2 = 1
for j in range(len(f)):
    if(abs(f[j] - 2*df) < e2):
        e2 = abs(f[j] - 2*df)
        pos2 = j


print(f[pos1])
print(f[pos2])

eff_4w = np.log10(np.abs(Azf[pos1])**2)-np.log10(np.abs(Azf[pos2])**2)
print(eff_4w)

plt.plot(t,np.abs(A0t)**2, 'k', linewidth = 0.8, label = 'Pulso inicial')
plt.plot(t,np.abs(Azt)**2, 'g', label = 'Pulso con XPM')
plt.plot(t,np.abs(Azt1)**2, 'r', label = 'Pulso sin XPM')
plt.title('Dominio del tiempo')
plt.axvline(8340*m.pi, ls = '--', color = 'g', linewidth = 0.5)
#plt.xlim(-120,120)
plt.xlabel('t[ps]')
plt.ylabel('|A(z, t)|^2 [u.a.]')
plt.legend()
plt.show()

plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(A0f)**2/0.001), 'k--', linewidth = 0.8, label = 'Pulso inicial') 
plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(Azf)**2/0.001), 'g', label = 'Pulso con XPM')
plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(Azf1)**2/0.001), 'r', label = 'Pulso sin XPM')
plt.title('Dominio de la frecuencia')
plt.xlabel('f[GHz]')
plt.ylabel('|A(z, f)|^2 [u.a.]')
#plt.xlim(-800,800)
plt.legend()
plt.show()


plt.plot(np.fft.fftshift(1000*f),10*np.log10(np.fft.fftshift(np.abs(A0f)**2)), 'k--', linewidth = 0.8, label = 'Pulso inicial') 
plt.plot(np.fft.fftshift(1000*f),10*np.log10(np.fft.fftshift(np.abs(Azf)**2)), 'g', label = 'Pulso con XPM')
plt.plot(np.fft.fftshift(1000*f),10*np.log10(np.fft.fftshift(np.abs(Azf1)**2)), 'r', label = 'Pulso sin XPM')
plt.title('Dominio de la frecuencia')
plt.xlabel('f[GHz]')
plt.ylabel('|A(z, f)|^2 [u.a en dB]')
#plt.xlim(-800,800)
#plt.ylim(-100)
plt.legend()
plt.show()
