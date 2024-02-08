from cProfile import label
from cmath import sqrt
from turtle import color
import scipy as sc
import matplotlib.pyplot as plt
import numpy as np
import lib
import math as m

#specs
b3 = 0.163 #ps^3/km

L = 31000 #largo de la fibra en km
d = 40#cant de pasos para el vector z
z = np.linspace(0, L, d)



n = 16384 #cantidad de muestras
dt = 0.05 #en picosegundos
f_s = 1/dt #en terahertz

t = np.linspace(-n*dt/2,n*dt/2,n-1)

f = np.linspace(-f_s/2,f_s/2,n-1)
f = sc.fft.fftshift(f)
f = np.roll(f,1)


t0 = 10 #en picosegundos
A0t = np.exp(-0.5*(t/t0)**2) #pulso inicial
A0f = sc.fft.fft(A0t) #su 

Ld = t0**3/b3

ancho0 = lib.ancho_RMS(t, np.abs(A0t)**2)

ancho = []
ancho_teo = []
for pos in z:
    Azf = A0f * np.exp((1/6)*1j*b3*pos*(2*np.pi*f)**3)
    Azt = sc.fft.ifft(Azf)
    ancho = np.append(ancho,lib.ancho_RMS(t, np.abs(Azt)**2)/ancho0)
    ancho_teo = np.append(ancho_teo,sqrt(1+0.5*(0.25*pos*b3/ancho0**3)**2))



plt.plot(t,np.abs(Azt)**2, 'r', label = 'Pulso en z = L')
plt.plot(t,np.abs(A0t)**2, 'k', label = 'Pulso en z = 0')
plt.title('Dominio del tiempo')
plt.xlabel('t[ps]')
plt.ylabel('|A(z, t)|^2')
plt.legend()
plt.show()

#shifteo todo para graficar los espectros correctamente
plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(A0f)**2), 'k', label = 'Pulso en z = 0') 
plt.plot(np.fft.fftshift(1000*f),np.fft.fftshift(np.abs(Azf)**2), 'r', label = 'Pulso en z = L')
plt.title('Dominio de la frecuencia')
plt.xlabel('f[GHz]')
plt.ylabel('|A(z, f)|^2')
plt.xlim(-100,100)
plt.legend()
plt.show()

plt.plot(z,ancho, 'bo', label = 'Experimental')
plt.plot(z,ancho_teo, 'r', label = 'Analítico')
plt.xlabel('z [km]')
plt.ylabel('Dispersión del pulso')
plt.axvline(5*Ld, ls = '--', color = 'k', label = "5 Ld'")
plt.legend()
plt.show()