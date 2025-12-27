# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 18:17:52 2024

@author: finnd
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import dblquad

#Define the real part of the integral ----------------------------------
def F2dreal(yp , xp , y , x , k, z):
    prefactor = (k*E_0)/(2*np.pi*z)
    kernel_r = np.cos(k/(2*z)*((x-xp)**2 +(y - yp)**2))
    return prefactor*kernel_r
    
# Define function for imaginary part of integral ---------------------------
def F2dimag( yp , xp , y , x , k ,z):
        prefactor = (k*E_0)/(2*np.pi*z)
        kernel_i = np.sin(k/(2*z)*((x-xp)**2 +(y - yp)**2))
        return prefactor*kernel_i
    

 
# Define the upper and lower limits for a square aperture-----------
apt_size = 1e-4
xp1 = -5e-5
xp2 = -xp1
yp1 = -xp1
yp2 = -yp1
#Define Physical Constants ---------------------------------
e_0 = 8.8e-12 #  epsilon 0 the permittivity of free space F/m
c = 3e8 # The speed of light m/s

z = 5e-2 # screen-aperture separation in m
k = (2*np.pi)/5e-9 # wave number in rad/m
E_0 = 1 # initial energy of wave in the aperture

obslim_x = 5e-2
obslim_y = 5e-2
num = 20

#y = -5e-5
#x=-5e-5




xvals = np.linspace(-obslim_x, obslim_x, num)
yvals = np.linspace(-obslim_y, obslim_y, num)
X,Y = np.meshgrid(xvals, yvals)
intensity = np.zeros((num,num))


for i in range(num):
    x = xvals[i]
    for j in range(num):
        
        y = yvals[j]
        realpart, realerror = dblquad(F2dreal, xp1, xp2, yp1, yp2, args = (y, x, k, z))
        imagpart, imagerror = dblquad(F2dimag, xp1, xp2, yp1, yp2, args = (y, x, k, z))     
        
        intensity[j, i] = e_0*c*(realpart**2 + imagpart**2)
        
    
       




#plt.imshow(intensity,  X, Y, norm = 'log', levels=100, cmap='plasma')
plt.imshow(intensity, norm = 'log', extent=[-obslim_x, obslim_x, -obslim_y, obslim_y], 
          cmap='plasma', origin='lower', aspect='auto')
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title('Rectangular diffraction')
plt.colorbar()
plt.xlim(-obslim_x, obslim_x)
plt.ylim(-obslim_y, obslim_y)

plt.show()




