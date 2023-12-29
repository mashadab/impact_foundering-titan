#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 08:39:07 2023

@author: afzal-admin
"""
import numpy as np
import matplotlib.pyplot as plt
#from skimage import color
#from skimage import io
from matplotlib import cm
from scipy.interpolate import interp2d, RegularGridInterpolator
import scipy.io as sio

image = plt.imread('./jpl_croppped.jpeg')

cmap_viridis = cm.get_cmap('viridis')


cmap_greys = cm.get_cmap('gray')


numbers = np.transpose([np.linspace(0,1,256)])
T_min = 90 #Minimum temperature [K]
T_max = 273.16 #Minimum temperature [K]

#getting the scales out
cmap_viridis_disc = (cmap_viridis(numbers))[:,0,0:3]
cmap_greys_disc   = (cmap_greys(numbers))[:,0,0]

T_disc   = T_min + (T_max - T_min) * (cmap_greys(numbers))[:,0,0]   #Temperature conversion [K]

#making one to one connection ~ Linear Algebra
conversion = np.linalg.inv(cmap_viridis_disc.T @ cmap_viridis_disc) @ cmap_viridis_disc.T  @ cmap_greys_disc

#converting image to grayscale 
image_rescale = image/255

phi = image_rescale[:,:,0]

phi[phi<0.9] = 0
phi[phi>=0.9] = 1
#image_rescale[image_rescale>1-1e-10] = #p.nan  #Getting rid of white lines

#summing the numbers
grayscale_image1 = 255*(conversion[0] * image_rescale[:,:,0] + conversion[1] * image_rescale[:,:,1] + conversion[2] * image_rescale[:,:,2])

grayscale_image = grayscale_image1[:,0:int(grayscale_image1.shape[1]/2)]

X = np.linspace(0,30e3,grayscale_image.shape[1])  #radius in meters
Y = np.linspace(20e3,0,grayscale_image.shape[0])   #height in meters

x, y = np.meshgrid(X, Y)

plt.imshow(grayscale_image, cmap='gray', vmin=0, vmax=255)



#Reducing the image size
compression_factor = 2
reduce_image = grayscale_image[::compression_factor,::2*compression_factor]
X_reduced = x[::compression_factor,::2*compression_factor]  #radius in meters
Y_reduced = y[::compression_factor,::2*compression_factor]  #height in meters

#melt fraction 
phi = reduce_image/255
phi[X_reduced >6.96e3] = 0.0 #Brute forcing 
phi[phi<0.9] = 0
phi[phi>=0.9] = 1

#Temperature [K]
Tdata = T_min + (T_max-T_min)*reduce_image/255


#Plotting

plt.figure()
plt.contourf(X_reduced,Y_reduced,reduce_image,cmap='gray')
plt.colorbar()

plt.figure()
plt.contourf(X_reduced,Y_reduced,Tdata,cmap='Reds', level=100,vmin=T_min,vmax=T_max)
plt.colorbar()

plt.figure()
plt.contourf(X_reduced,Y_reduced,phi,cmap='Blues', level=100,vmin=0,vmax=1)
plt.colorbar()



matlab_data = dict(phidata=phi,Xdata=X_reduced,Ydata=Y_reduced,Tdata=Tdata)
sio.savemat('data_dump.mat', matlab_data)
