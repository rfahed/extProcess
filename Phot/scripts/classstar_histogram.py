#!/usr/bin/env python

from __future__ import division
from Phot import catalog
import numpy as np
import matplotlib.pyplot as P

filter ="g"

incat="input_{}_matchtag.txt".format(filter)
outcat="output_{}_seeing_1.txt".format(filter)
mergedcat="merged_{}.txt".format(filter)
realcat="real_{}.txt".format(filter)

inc=catalog.read(incat)
outc=catalog.read(outcat)
merc=catalog.read(mergedcat)
realc=catalog.read(realcat)

# filter bad detection (close to saturation or border)
merc=merc[merc['FLAGS'] == 0]
outc=outc[outc['FLAGS'] == 0]
realc=realc[realc['FLAGS'] == 0]

stars = merc[merc['OBJTYPE']==0]
gals = merc[merc['OBJTYPE']>0]

if False:
    minb=0
    maxb=1
    binw=0.1
    ylim = [0,25000]
    bins=np.linspace(minb,maxb,int((maxb-minb)/binw))

    P.subplot(2, 1, 1)
    n_in,bins,patches_in=catalog.histogramcol(stars,"CLASS_STAR",label="stars", bins=bins)
    P.legend()
    P.ylim(ylim)
    
    P.subplot(2, 1, 2)
    n_out,bins,patches_out=catalog.histogramcol(gals,"CLASS_STAR", label="galaxies", bins=bins)
    P.legend()
    P.ylim(ylim)
    P.show()

if False :
    minb=10
    maxb=30
    binw=1

    bins=np.linspace(minb,maxb,int((maxb-minb)/binw)+1)
    sexstars = outc[outc['CLASS_STAR'] > 0.8]
    sexgals = outc[outc['CLASS_STAR'] < 0.2] 
    unknowns = outc[(outc['CLASS_STAR'] > 0.2) & (outc['CLASS_STAR'] < 0.8)] 
    
    n,bins=np.histogram(outc["MAG_AUTO"], bins=bins)
    n_stars,bins=np.histogram(sexstars["MAG_AUTO"], bins=bins)
    n_gals,bins=np.histogram(sexgals["MAG_AUTO"], bins=bins)
    n_unk,bins=np.histogram(unknowns["MAG_AUTO"], bins=bins)
    
    mags=0.5*(bins[:-1]+bins[1:])
    
    P.figure()
    P.plot(mags,n_stars/n,label="fraction of stars ({} band)".format(filter))
    P.plot(mags,n_gals/n,label="fraction of gals ({} band)".format(filter))
    P.plot(mags,n_unk/n,label="fraction unidentified objects ({} band)".format(filter))
    P.legend(loc=0)
    P.xlabel('Magnitude')
    P.ylabel('Fraction')
    P.ylim([-0.1,1.1])
    P.xlim([10,26])
    P.tight_layout()
    P.savefig("class_fraction_sim_{}.png".format(filter))
#    P.show()

if False :
    minb=10
    maxb=30
    binw=1

    bins=np.linspace(minb,maxb,int((maxb-minb)/binw)+1)
    sexstars = realc[realc['CLASS_STAR'] > 0.8]
    sexgals = realc[realc['CLASS_STAR'] < 0.2] 
    unknowns = realc[(realc['CLASS_STAR'] > 0.2) & (realc['CLASS_STAR'] < 0.8)] 
    
    n,bins=np.histogram(realc["MAG_AUTO"], bins=bins)
    n_stars,bins=np.histogram(sexstars["MAG_AUTO"], bins=bins)
    n_gals,bins=np.histogram(sexgals["MAG_AUTO"], bins=bins)
    n_unk,bins=np.histogram(unknowns["MAG_AUTO"], bins=bins)
    
    mags=0.5*(bins[:-1]+bins[1:])

    P.figure()
    P.plot(mags,n_stars/n,label="fraction of stars ({} band)".format(filter))
    P.plot(mags,n_gals/n,label="fraction of gals ({} band)".format(filter))
    P.plot(mags,n_unk/n,label="fraction unidentified objects ({} band)".format(filter))
    P.legend(loc=0)
    P.xlabel('Magnitude')
    P.ylabel('Fraction')
    P.ylim([-0.1,1.1])
    P.xlim([10,26])
    P.tight_layout()
    P.savefig("class_fraction_real_{}.png".format(filter))

if False :
    P.figure()
    catalog.scattercols(realc,"MAG_AUTO","CLASS_STAR")
    P.xlim([10,30])
    P.grid()
    P.savefig("classstar_mag_real_{}.png".format(filter))

if False :
    P.figure()
    #catalog.scattercols(stars,"MAG_AUTO","CLASS_STAR",color="r",label="stars")
    #catalog.scattercols(gals,"MAG_AUTO","CLASS_STAR",color="b",label="gals")
    catalog.scattercols(outc,"MAG_AUTO","CLASS_STAR")
    P.grid()
    P.xlim([10,30])
    P.legend()
    P.savefig("classstar_mag_real_{}.png".format(filter))
    P.show()
    
if False :
    minb=0
    maxb=5
    binw=0.1

    bins=np.linspace(minb,maxb,int((maxb-minb)/binw)+1)
    
    P.figure()
    stars["FWHM_WORLD"] = stars["FWHM_WORLD"]*3600.
    catalog.histogramcol(stars,"FWHM_WORLD",bins=bins,label="stars in sim image ({} band)".format(filter))
    P.grid()
    P.show()

if True :
    minb=0
    maxb=5
    binw=0.1

    bins=np.linspace(minb,maxb,int((maxb-minb)/binw)+1)
    
    P.figure()
    realc["FWHM_WORLD"] = realc["FWHM_WORLD"]*3600.
    sexstars = realc[realc['CLASS_STAR'] > 0.8]
    sexgals = realc[realc['CLASS_STAR'] < 0.2]
    unknowns = realc[(realc['CLASS_STAR'] > 0.2) & (realc['CLASS_STAR'] < 0.8)]
    catalog.histogramcol(sexstars,"FWHM_WORLD",bins=bins,label="stars in real image ({} band)".format(filter))
    catalog.histogramcol(sexgals,"FWHM_WORLD",bins=bins,label="galaxies in real image ({} band)".format(filter))
    catalog.histogramcol(unknowns,"FWHM_WORLD",bins=bins,label="undefined obj in real image ({} band)".format(filter))
    P.grid()
    P.show()
