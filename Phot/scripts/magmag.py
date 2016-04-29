#!/usr/bin/env python

from astropy.io import ascii
from Phot import catalog

catalogs = ['mag_DESg_gal.txt','mag_DESr_gal.txt','mag_DESi_gal.txt','mag_DESz_gal.txt','mag_DESY_gal.txt']

cats = []
for cat in catalogs :
    c=ascii.read(cat)
    catalog
    

