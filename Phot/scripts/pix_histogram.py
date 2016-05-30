#!/usr/bin/env python

from Phot import catalog
from Phot import image
import numpy as np
from astropy.io import fits 
from matplotlib import pylab as p
from scipy.optimize import curve_fit

filter="g"
sim_image="DES{}_sim.fits".format(filter)
real_image="DES_{}.fits".format(filter)
ext=1
im_sim=fits.open(sim_image)
pix_sim=np.ravel(im_sim[ext].data)

im_real=fits.open(real_image)
pix_real=np.ravel(im_real[ext].data)

nbins=1301
minb=-30
bins=np.linspace(minb,minb+nbins-1,nbins)
n_sim, bins, patches = p.hist(pix_sim, bins=bins, alpha=0.5,label="sim image (band {})".format(filter))
n_real, bins, patches = p.hist(pix_real, bins=bins, alpha=0.5,label="real image (band {})".format(filter))
p.xlabel('pixel value (ADU)',size=20)
p.ylabel('counts',size=20)


if True :
    bins = (bins[:-1] + bins[1:])/2

    i_s=np.where(n_sim==max(n_sim))
    i_s=i_s[0][0]
    p0_s=[5000,bins[i_s],10]

    p_s, var_matrix = curve_fit(image.gauss, bins, n_sim, p0=p0_s)

    i_r=np.where(n_real==max(n_real))
    i_r=i_r[0][0]
    p0_r=[5000,bins[i_r],10]

    p_r, var_matrix = curve_fit(image.gauss, bins, n_real, p0=p0_r)

    gauss_s = image.gauss(bins,*p_s)
    gauss_r = image.gauss(bins,*p_r)

    p.plot(bins, gauss_s, '--',label='sigma = {}'.format(p_s[2]),linewidth=2)
    p.plot(bins, gauss_r, '--',label='sigma = {}'.format(p_r[2]),linewidth=2)

p.grid()
p.legend()
p.show()

