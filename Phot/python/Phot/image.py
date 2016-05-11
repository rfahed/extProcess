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
import numpy as np


def sex(imname, zeropoint=0):

    catalog_name = utils.rm_extension(imname)+".txt"
    cmd = ["sex",imname,"-c",utils.getAuxPathFile("default.sex")]
    cmd += ["-MAG_ZEROPOINT",str(zeropoint)]
    cmd += ["-CATALOG_NAME",catalog_name]
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


