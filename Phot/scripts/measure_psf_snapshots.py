#!/usr/bin/env python

from astropy.io import fits
import sys
from Phot import image
import numpy as np
from matplotlib import pylab as p

def measure_psf_snapshots(snapcat):
    snapcat=fits.open(snapcat)
    #pixscale=image.get_pixscale(snapcat[0].header)
    pixscale=0.19*0.263
    shape=snapcat[0].data.shape
    vignets=snapcat[0].data.reshape(shape[0]*shape[1],shape[2],shape[3])

    psfs=image.measure_psfs(vignets, pixscale=pixscale, plot=True)
    return np.median(psfs)
    
if __name__ == "__main__" :
    psfs = []
    for psffile in sys.argv[1:]:
        psf = measure_psf_snapshots(psffile)
        psfs.append(psf)
        print psf
    
    if len(psfs) > 10:    
        p.hist(psfs,normed=True)
        p.xlabel('psf size')
        p.tight_layout()
        p.savefig("psfsize.png")