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

def get_wcs(header):
    KEYS_TO_DEL = ["PV{}_{}".format(i,j) for i in range(11) for j in range(11) ]
    for k in KEYS_TO_DEL:
        if k in header :
            del header[k]

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
    if show :
        P.show()
   

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def sex(imname, zeropoint=0,outputcat=None):
    with open(utils.getAuxPathFile("default.sex.template"),'r') as f :
        conftemp = Template(f.read())
    conf = conftemp.substitute({'zeropoint':str(zeropoint)})
    with open("default.sex",'w') as f :
	f.write(conf)
    if outputcat is None :
    	outputcat = utils.rm_extension(imname)+".txt"

    cmd = ["sex",imname,"-c","default.sex"]
    cmd += ["-CATALOG_NAME",outputcat]
    cmd += ["-PARAMETERS_NAME",utils.getAuxPathFile("default.param")]
    cmd += ["-FILTER_NAME",utils.getAuxPathFile("gauss_5.0_9x9.conv")]
    cmd += ["-STARNNW_NAME",utils.getAuxPathFile("default.nnw")]

    p = subprocess.call(cmd)
        
    if p!=0 :
        sys.exit("SExtractor failed... Exiting.")

def get_zeropoint(imname,apply_exptime=False,zerokey="MAGZERO"):
    im = fits.open(imname)
    try :
	zeropoint = im[0].header[zerokey]
    except :
	zeropoint = im[1].header[zerokey]

    if apply_exptime :
	zeropoint += 2.5*np.log10(im[0].header['EXPTIME'])

    logger.info("Mag zeropoint : {}".format(zeropoint))
    return zeropoint


