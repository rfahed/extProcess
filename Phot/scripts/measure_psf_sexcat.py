#!/usr/bin/env python

from astropy.io import fits
import sys
from Phot import image
import numpy as np



def measure_psf_sexcat(sexcat):
    sexcat=fits.open(sexcat)
    h=fits.Header.fromstring("\n".join(sexcat[1].data[0][0]), sep="\n")
    pixscale=image.get_pixscale(h)
    vignets=sexcat[2].data["VIGNET"]

    psfs=image.measure_psfs(vignets, pixscale=pixscale, plot=False, mask_value=-1e+30)
    return np.median(psfs)
    
if __name__ == "__main__" :
    psfs = []
    for psffile in sys.argv[1:]:
        psf = measure_psf_sexcat(psffile)
        psfs.append(psf)
        print psf
    
    if len(psfs) > 10:    
        p.hist(psfs,normed=True)
        p.xlabel('psf size')
        p.tight_layout()
        p.savefig("psfsize.png")