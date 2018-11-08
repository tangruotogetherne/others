#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 10:44:22 2018

@author: tac
"""

# Surface Brightness Limit

# <------ Step 1. Get the initial background value ------>


import numpy as np
import astropy.io.fits as fits
from photutils import CircularAperture, aperture_photometry      # Circular aperture is more exact 
                                                                 # than rectangular shape on photutils
from astropy.visualization import ImageNormalize, ZScaleInterval
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats

matplotlib.rc('font', family='serif', size=15)
matplotlib.rc('xtick', direction='in', labelsize=12)
matplotlib.rc('ytick', direction='in', labelsize=12)

pixel = 0.186
x_arcsec = 10
r = np.sqrt((x_arcsec/pixel)**2/(np.pi))
photoz = 30.

hdu_seg = fits.open('/home/tac/Documents/SBlimit/new/D1.R.seg03.fits')               # Segmentation image
hdu_bkg = fits.open('/home/tac/Documents/SBlimit/new/D1.R.bkg03_MaskedObjects.fits') # Background image
seg_data = hdu_seg[0].data
bkg_data = hdu_bkg[0].data

# Define the positions of apertures with 30000 2D random data
# Avoiding the positions of none value
# The total pixel size is pixels_x = 20347, pixels_y = 21404

n = 30000
positions = np.eye(n,2)
x_min = 360     
x_max = 19720
y_min = 620
y_max = 20594
x = np.random.randint(x_min+r, x_max-r+1, size=(n,1))
y = np.random.randint(y_min+r, y_max-r+1, size=(n,1))
positions[:,0] = x[:,0]
positions[:,1] = y[:,0]

# Aperture photometry for determining flux value
# Record positions with theirselves counts as flux value into table_new.txt

flux = np.zeros((n,1))
pos = np.zeros((n,2))
table = np.zeros((n,3))
i = 0
while i < n:
    apertures = CircularAperture(positions[i], r)
    seg_table = aperture_photometry(seg_data, apertures)
    counts = seg_table['aperture_sum']
    if counts > 0:
        flux[i] = -10**8
        pos[i] = positions[i]
    else:
        bkg_table = aperture_photometry(bkg_data, apertures)
        flux[i] = bkg_table['aperture_sum']
        pos[i] = positions[i]
    i+=1
table[:,0] = pos[:,0]
table[:,1] = pos[:,1]
table[:,2] = flux[:,0]
np.savetxt('table_new.txt', table)

# Examine the positions got from random data
# Whether the distribution is uniform or not 

data = np.loadtxt('table_new.txt')
x = data[:,0]
xgrid = np.arange(360, 19720, 1000)
fig = plt.figure(figsize=[8,8])
ax = fig.add_subplot(111)
ax.hist(x, xgrid, edgecolor='k', alpha=0.4)
ax.set_xlabel('Positions')
ax.set_ylabel('Counts')
plt.savefig('pos_wave')

# <------ Step 2. Choose the background apertures and extract more accurate value ------>


# Filtering the bad value
# Calculate a initial standard deviation and SBlimit
# Print the filtered results and the length of each result

arr = []
for value in flux:
    if value > -10**8:
        arr.append(value)
flux_0 = np.array(arr)
np.savetxt('flux0.txt', flux_0)
print('length of flux0', len(flux_0))
sdev_0 = np.std(flux_0)
print('sigma of flux0', sdev_0)
M_lim = -2.5*np.log10(sdev_0/(np.pi*r**2*pixel**2)) + photoz
print('M_lim_0 =', M_lim, '\n')
f = open('SBlimit_print_new.txt','w+')
f.write('\nlength of flux0 ')
f.write(str(len(flux_0)))
f.write('\nsigma of flux0 ')
f.write(str(sdev_0))
f.write('\nM_lim_0 = ')
f.write(str(M_lim))
f.write('\n')

# Clipping with 3-sigma clipping

# The first clipping

arr = []
sig = 3*sdev_0
i = 0
while i < len(flux_0):
    bkg = np.abs(flux_0[i] - np.mean(flux_0))
    if bkg < sig:
        arr.append(flux_0[i])
    i+=1
flux_1 = np.array(arr)
np.savetxt('flux1.txt', flux_1)
print('length of flux1', len(flux_1))
sdev_1 = np.std(flux_1)
print('sigma of flux1', sdev_1)
M_lim = -2.5*np.log10(sdev_1/(np.pi*r**2*pixel**2)) + photoz
print('M_lim_1 =', M_lim, '\n')
f.write('\nlength of flux1 ')
f.write(str(len(flux_1)))
f.write('\nsigma of flux1 ')
f.write(str(sdev_1))
f.write('\nM_lim_1 = ')
f.write(str(M_lim))
f.write('\n')

# The second clipping

arr = []
sig = 3*sdev_1
i = 0
while i < len(flux_1):
    bkg = np.abs(flux_1[i] - np.mean(flux_1))
    if bkg < sig:
        arr.append(flux_1[i])
    i+=1
flux_2 = np.array(arr)
np.savetxt('flux2.txt', flux_2)
print('length of flux2', len(flux_2))
sdev_2 = np.std(flux_2)
print('sigma of flux2', sdev_2)
M_lim = -2.5*np.log10(sdev_2/(np.pi*r**2*pixel**2)) + photoz
print('M_lim_2 =', M_lim, '\n')
f.write('\nlength of flux2 ')
f.write(str(len(flux_2)))
f.write('\nsigma of flux2 ')
f.write(str(sdev_2))
f.write('\nM_lim_2 = ')
f.write(str(M_lim))
f.write('\n')

# The third clipping

arr = []
sig = 3*sdev_2
i = 0
while i < len(flux_2):
    bkg = np.abs(flux_2[i] - np.mean(flux_2))
    if bkg < sig:
        arr.append(flux_2[i])
    i+=1
flux_3 = np.array(arr)
np.savetxt('flux3.txt', flux_3)
print('length of flux3', len(flux_3))
sdev_3 = np.std(flux_3)
print('sigma of flux3', sdev_3)
M_lim = -2.5*np.log10(sdev_3/(np.pi*r**2*pixel**2)) + photoz
print('M_lim_3 =', M_lim, '\n')
f.write('\nlength of flux3 ')
f.write(str(len(flux_3)))
f.write('\nsigma of flux3 ')
f.write(str(sdev_3))
f.write('\nM_lim_3 = ')
f.write(str(M_lim))
f.write('\n')

# The fourth clipping

arr = []
sig = 3*sdev_3
i = 0
while i <len(flux_3):
    bkg = np.abs(flux_3[i] - np.mean(flux_3))
    if bkg < sig:
        arr.append(flux_3[i])
    i+=1
flux_4 = np.array(arr)
np.savetxt('flux4.txt', flux_4)
print('length of flux4', len(flux_4))
sdev_4 = np.std(flux_4)
print('sigma of flux4', sdev_4)
M_lim = -2.5*np.log10(sdev_4/(np.pi*r**2*pixel**2)) + photoz
print('M_lim_4 =', M_lim, '\n')
f.write('\nlength of flux4 ')
f.write(str(len(flux_4)))
f.write('\nsigma of flux4 ')
f.write(str(sdev_4))
f.write('\nM_lim_4 = ')
f.write(str(M_lim))
f.write('\n')

# The fifth clipping

arr = []
sig = 3*sdev_4
i = 0
while i < len(flux_4):
    bkg = np.abs(flux_4[i] - np.mean(flux_4))
    if bkg < sig:
        arr.append(flux_4[i])
    i+=1
flux_5 = np.array(arr)
np.savetxt('flux5.txt', flux_5)
print('length of flux5', len(flux_5))
sdev_5 = np.std(flux_5)
print('sigma of flux5', sdev_5)
M_lim = -2.5*np.log10(sdev_5/(np.pi*r**2*pixel**2)) + photoz
M = M_lim
print('M_lim_5 =', M_lim, '\n')
f.write('\nlength of flux5 ')
f.write(str(len(flux_5)))
f.write('\nsigma of flux5 ')
f.write(str(sdev_5))
f.write('\nM_lim_5 = ')
f.write(str(M_lim))
f.write('\n')

# The sixth clipping

arr = []
sig = 3*sdev_5
i = 0
while i < len(flux_5):
    bkg = np.abs(flux_5[i] - np.mean(flux_5))
    if bkg < sig:
        arr.append(flux_5[i])
    i+=1
flux_6 = np.array(arr)
np.savetxt('flux6.txt', flux_6)
print('length of flux6', len(flux_6))
sdev_6 = np.std(flux_6)
print('sigma of flux6', sdev_6)
M_lim = -2.5*np.log10(sdev_6/(np.pi*r**2*pixel**2)) + photoz
print('M_lim_6 =', M_lim, '\n')
f.write('\nlength of flux6 ')
f.write(str(len(flux_6)))
f.write('\nsigma of flux6 ')
f.write(str(sdev_6))
f.write('\nM_lim_6 = ')
f.write(str(M_lim))
f.write('\n')

# <------ Step 3. Show the apertures of last value at bkg image ------>


# Select flux_5 as the last result whose SBlimit and sigma all became an constant respectively
# Get every positions of value from the last flux array
# Paint the last apertures at bkg image

x = []
y = []
i = 0
while i < len(table):
    j = 0
    while j < len(flux_5):
        if table[:,2][i] == flux_5[j]:
            x.append(table[:,0][i])
            y.append(table[:,1][i])
        j+=1
    i+=1
l = len(x)
positions_filter = np.zeros((l,2))
positions_filter[:,0] = x
positions_filter[:,1] = y
np.savetxt('pos_fil_new.txt', positions_filter)

print('The number of final value :', len(positions_filter))
f.write('\nThe number of final value : ')
f.write(str(len(positions_filter)))
f.write('\n')

norm = ImageNormalize(bkg_data, ZScaleInterval())
fig = plt.figure(figsize=[8,8])
apertures = CircularAperture(positions_filter, r=r)
apertures.plot(color='r', lw=1.5, alpha=1.)
plt.imshow(bkg_data, norm=norm, origin='lower', cmap='gray_r') # Accroding to this distribution of apertures
                                                               # it could determine the region of local image
plt.savefig('all_aper_new')

# (optional) Show the apertures' positions distribution with histogram

pos_x = positions_filter[:,0]
fig = plt.figure(figsize=[8,8])
ax = fig.add_subplot(111)
xgrid = np.arange(360, 19720, 1000)
ax.hist(pos_x, xgrid, edgecolor='k', alpha=0.4)
ax.set_title('Distribution of filtered positions')
ax.set_ylabel('Counts')
ax.set_xlabel('Positions')
plt.savefig('pos_fil_wave')

# <------ Step 4. Calculate background variation ------>


# Show the distribution of background with histogram
# Using probability density function of norm for fitting

xgrid = np.arange(-800,200,30) 
xcenter = (xgrid[1:]+xgrid[:-1])/2 
hx, xedge = np.histogram(flux_5, xgrid)

mean = np.mean(flux_5) 
sigma = np.std(flux_5) 
norm_pdf = stats.norm.pdf(xcenter, loc = mean, scale = sigma)  # the probability density function of
                                                               # normal distribution

print('flux_5 :', mean, sigma) 
f.write('\nflux_5 ')
f.write(str(mean))
f.write(' ')
f.write(str(sigma))
f.write('\n')

fig = plt.figure(figsize=[10,8])
ax = fig.add_subplot(111)
ax.plot(xcenter, norm_pdf/np.sum(norm_pdf)*np.sum(hx), 'r-', markersize=20, label='normal distribution') 
ax.hist(flux_5, xgrid, edgecolor='k', facecolor='w', label='background')
legend = plt.legend(fontsize='smaller')
plt.savefig('fitting_image_new')

# <------ Step 5. Comparison the different type of aperture ------>


# Draw the apertures containing objects and apetures just for bkg

x = []
y = []
i = 0
while i < len(table):
    if table[:,2][i] == -10**8:
        x.append(table[:,0][i])
        y.append(table[:,1][i])
    i+=1
m = len(x)
pos_obj = np.zeros((m,2))
pos_obj[:,0] = x
pos_obj[:,1] = y
pos_cut1 = np.zeros((l,2))
pos_cut2 = np.zeros((m,2))
pos_cut1[:,0] = positions_filter[:,0]-6000  # The coordinate of apertures must follow the real frame of image
pos_cut1[:,1] = positions_filter[:,1]-2700
pos_cut2[:,0] = pos_obj[:,0]-6000
pos_cut2[:,1] = pos_obj[:,1]-2700

fig, ax = plt.subplots(figsize=[8,8])
apertures1 = CircularAperture(pos_cut1, r=r)
apertures2 = CircularAperture(pos_cut2, r=r)
apertures1.plot(color='r', lw=1.5, alpha=1.)
apertures2.plot(color='b', lw=1.5, alpha=1.)
plt.imshow(bkg_data[2700:3200, 6000:6500], norm=norm, origin='lower', cmap='gray_r')
plt.savefig('comparison_image_new')
print('the number of red apertures is ', l , '\nthe number of blue apertures is ', m)
print('the surface brightness limit is ', M)
f.write('\nthe number of red apertures is ')
f.write(str(l))
f.write('\nthe number of blue apertures is ')
f.write(str(m))
f.write('\n')
f.write('\nthe surface brightness limit is ')
f.write(str(M))
f.close()