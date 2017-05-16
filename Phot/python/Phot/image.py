"""
@file python/Phot/photometry.py
@date 03/16/16
@author user
"""
import ElementsKernel.Logging as log
logger = log.getLogger('image')
import subprocess, os, sys
from . import utils
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
import numpy as np
from string import Template
from matplotlib import pylab as P
import json

def get_wcs(header):
    #KEYS_TO_DEL = ["PV{}_{}".format(i,j) for i in range(11) for j in range(11) ]
    #for k in KEYS_TO_DEL:
    #    if k in header :
    #        del header[k]

    wcs = WCS(header)
    return wcs
    
def get_footprint(header):
    wcs = get_wcs(header)
    footprint = wcs.calc_footprint(header)
    return footprint

def in_rect(rect,x,y):
    # Corners in ax,ay,bx,by,dx,dy
    # Point in x, y
    ax=rect[0][0]
    ay=rect[0][1]
    bx=rect[1][0]
    by=rect[1][1]
    dx=rect[3][0]
    dy=rect[3][1]

    bax = bx - ax
    bay = by - ay
    dax = dx - ax
    day = dy - ay

    if ((x - ax) * bax + (y - ay) * bay < 0.0):
        return False
    if ((x - bx) * bax + (y - by) * bay > 0.0):
        return False
    if ((x - ax) * dax + (y - ay) * day < 0.0): 
        return False
    if ((x - dx) * dax + (y - dy) * day > 0.0): 
        return False

    return True

def plotfits(imname,ext=1,show=False,**kwargs):
    hdu = fits.open(imname)[ext]
    wcs = get_wcs(hdu.header)
    fig = P.figure()
    fig.add_subplot(111, projection=wcs)
    P.imshow(hdu.data, origin='lower', norm=LogNorm(), cmap='gray_r',**kwargs)
    P.colorbar()
    P.xlabel('RA')
    P.ylabel('Dec')

    return wcs

def plotcat(cat,wcs,show=False,poskeys=['X_WORLD','Y_WORLD'],**kwargs):
    objects_pos=wcs.all_world2pix(cat[poskeys[0]], cat[poskeys[1]], 0)
    P.scatter(objects_pos[0],objects_pos[1],**kwargs)
    if show :
        P.show()
    return objects_pos

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def smag2pix(smag,pixscale,zp,exptime=1.):
    pixvalue=exptime*10.**((zp-smag)/2.5)*pixscale**2
    return pixvalue

def sigma_background(counts_bg, readout_noise, gain):
    return np.sqrt(1/gain*(counts_bg + readout_noise))

def read_instrument(instrument_file):
    with open(instrument_file) as f :
        instrument = json.load(f)
    instrument["CCDS_LAYOUT"] = {int(k):v for (k,v) in instrument["CCDS_LAYOUT"].items()}
    return instrument

def updatewcs(im,instrument):
    for i,ext in enumerate(im[1:],1):
        wcs=instrument['CCDS_LAYOUT'][i]["WCS"]
        for k,v in wcs.items():
            ext.header[k]=v

def sex(imname, zeropoint=0,outputcat=None):
    with open(utils.getAuxPathFile("default.sex.template"),'r') as f :
        conftemp = Template(f.read())
    conf = conftemp.substitute({'zeropoint':str(zeropoint)})
    with open("default.sex",'w') as f :
        f.write(conf)
    if outputcat is None :
        outputcat = utils.rm_extension(imname)+".cat"

    cmd = ["sex",imname,"-c","default.sex"]
    cmd += ["-CATALOG_NAME",outputcat]
    cmd += ["-PARAMETERS_NAME",utils.getAuxPathFile("default.param")]
    cmd += ["-FILTER_NAME",utils.getAuxPathFile("gauss_5.0_9x9.conv")]
    cmd += ["-STARNNW_NAME",utils.getAuxPathFile("default.nnw")]

    p = subprocess.call(cmd)
        
    if p!=0 :
        sys.exit("SExtractor failed... Exiting.")

def get_zeropoints(imname,apply_exptime=False,zerokey="SIMMAGZP"):
    im = fits.open(imname)
    zeropoints=[]
    for ext in im[1:]:
        z=ext.header[zerokey]
        if apply_exptime :
            z += 2.5*np.log10(im[0].header['EXPTIME'])
        zeropoints.append(z)

    logger.info("Mag zeropoints : {}".format(zeropoints))
    return zeropoints


