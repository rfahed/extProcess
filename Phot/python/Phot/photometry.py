"""
@file python/Phot/photometry.py
@date 03/16/16
@author user
"""
import subprocess, os, sys
from . import utils    

def sex(imname):
    
    catalog_name = utils.rm_extension(imname)+".txt"
    cmd = ["sex",imname,"-c",utils.getAuxPathFile("default.sex")]
    cmd += ["-CATALOG_NAME",catalog_name]
    cmd += ["-PARAMETERS_NAME",utils.getAuxPathFile("default.param")]
    cmd += ["-FILTER_NAME",utils.getAuxPathFile("gauss_5.0_9x9.conv")]
    cmd += ["-STARNNW_NAME",utils.getAuxPathFile("default.nnw")]

    p = subprocess.call(cmd)
        
    if p!=0 :
        sys.exit("SExtractor failed... Exiting.")
