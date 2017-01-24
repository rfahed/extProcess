"""
@file python/Phot/catalog.py
@date 03/18/16
@author user
"""

import ElementsKernel.Logging as log
logger = log.getLogger('catalog')
from . import utils
from astropy.coordinates import ICRS, SkyCoord
from astropy import table
from astropy import units as u
from astropy.io import ascii,fits
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as P
import sys


NEWLINE = '\n'
MAG_ZEROPOINT = 30.5862157412

## Extension number of output catalog corresponds to CCD number
list_ext_num = [1, 2, 3, 4, 17, 18, 19, 20, 5, 6, 7, 8, 21, 22, 23, 24, 9, 10, 11,
            12, 25, 26, 27, 28, 13, 14, 15, 16, 29, 30, 31, 32]

def mergecats(catalogs=None,delta=1e-4,filters=None,poskeys=['X_WORLD','Y_WORLD'],stack=True,mcats=None):
    ncat = len(catalogs)
    if ncat > 2 :
        mcats=None
    if filters and stack :
        _poskeys = []
        for i in range(ncat) :
            catalogs[i] = tag_catalog(catalogs[i], filters[i])
            _poskeys.append([x+'_'+filters[i] for x in poskeys])
      
    mergedcat = catalogs[0]

    for i in range(1,ncat):
        if filters :
            mergedcat=mergecat(mergedcat,catalogs[i],poskeys1=_poskeys[0],poskeys2=_poskeys[i],delta=delta,stack=stack,mcats=mcats)
        else :
            mergedcat=mergecat(mergedcat,catalogs[i], poskeys1=poskeys, poskeys2=poskeys, delta=delta,stack=stack,mcats=mcats)
            
    return mergedcat
    
def plot_vignet(cat,i,show=True):
    P.imshow(cat['VIGNET'][i],cmap='gray_r',interpolation='none')
    P.colorbar()
    if show :
        P.show()

def add_calc_column(catalog,field1,field2,func,outputfield=None):
    if outputfield is None:
        outputfield = "calc_{}_{}".format(field1,field2)
    diff = table.Column(name=outputfield,data=func(catalog[field1],catalog[field2]))
    catalog.add_column(diff)

def add_diff_column(catalog,field1,field2,outputfield=None):
    if outputfield is None:
        outputfield = "diff_{}_{}".format(field1,field2)
    diff = table.Column(name=outputfield,data=catalog[field1]-catalog[field2])
    catalog.add_column(diff)

def add_sum_column(catalog,field1,field2,outputfield=None):
    if outputfield is None:
        outputfield = "sum_{}_{}".format(field1,field2)
    diff = table.Column(name=outputfield,data=catalog[field1]+catalog[field2])
    catalog.add_column(diff)

def add_frac_column(catalog,field1,field2,outputfield=None):
    if outputfield is None:
        outputfield = "frac_{}_{}".format(field1,field2)
    diff = table.Column(name=outputfield,data=catalog[field1]/catalog[field2])
    catalog.add_column(diff)

def add_prod_column(catalog,field1,field2,outputfield=None):
    if outputfield is None:
        outputfield = "prod_{}_{}".format(field1,field2)
    diff = table.Column(name=outputfield,data=catalog[field1]*catalog[field2])
    catalog.add_column(diff)

def plotcols(catalog,field1,field2,title=None,xlab=None,ylab=None,show=False,p=P,**kwargs):      
    p.plot(catalog[field1],catalog[field2],**kwargs)
    p.title(title)
    if xlab is None :
        xlab = field1
    if ylab is None :
        ylab = field2
    p.xlabel(xlab)
    p.ylabel(ylab)
    p.grid()
    if show:
        p.show()
        
def scattercols(catalog,field1,field2,title=None,xlab=None,ylab=None,log=False,show=False,p=P,**kwargs):    
    p.scatter(catalog[field1],catalog[field2],marker='+',**kwargs)
    if title :    
        p.title(title)
    if xlab is None :
        xlab = field1
    if ylab is None :
        ylab = field2
    if log:
        ax=p.gca()
        ax.set_yscale('log')

    p.xlabel(xlab)
    p.ylabel(ylab)
    p.grid()
    if show:
        p.show()
    
def histogramcol(catalog,field,xlab=None,ylab="Counts",show=False,p=P,**kwargs):
    if xlab is None :
        xlab = field
    n, bins, patches = p.hist(np.array(catalog[field]), alpha=0.5,**kwargs)
    p.xlabel(xlab)
    p.ylabel(ylab)
    p.legend(loc=0)
    p.grid()
    if show:
        p.show()

    return n, bins, patches    

def tag_catalog(catalog, tag):
    for cname in catalog.colnames :
        catalog.rename_column(cname,cname+'_'+tag)
    return catalog

def matchcats(p1,p2,delta=1e-4):
    #index,dist2d,dist3d = p1.match_to_catalog_sky(p2)
    #imatch = np.where(dist2d < delta*u.degree)
    imatch2, imatch1, d2d, d3d = p1.search_around_sky(p2, delta*u.deg)
    #return (imatch, index[imatch])
    return (imatch1,imatch2,d2d)

def mergecat(catalog1,catalog2, delta=1e-4, poskeys1=['X_WORLD','Y_WORLD'], poskeys2=['X_WORLD','Y_WORLD'], mcats=['cat1_w_matchtags.cat','cat2_w_matchtags.cat'], stack=True):    
    p1=SkyCoord(catalog1[poskeys1[0]],catalog1[poskeys1[1]],frame='icrs',unit="deg")
    p2=SkyCoord(catalog2[poskeys2[0]],catalog2[poskeys2[1]],frame='icrs',unit="deg")
    
    imatch1, imatch2, d = matchcats(p1,p2,delta=delta)
    distance=table.Column(name='DISTANCE',data=d) * 3600.0
    m1=np.zeros(len(catalog1),dtype=int)
    m1[imatch1]=1
    m2=np.zeros(len(catalog2),dtype=int)
    m2[imatch2]=1
    matched1=table.Column(name='MATCHED',data=m1)
    matched2=table.Column(name='MATCHED',data=m2)

    if mcats :
        catalog1.add_column(matched1)
        catalog2.add_column(matched2)
        #with open(mcats[0], 'w') as f :
            #ascii.write(catalog1, f,Writer=ascii.CommentedHeader)
        writefits(catalog1,mcats[0])
        #with open(mcats[1], 'w') as f :
            #ascii.write(catalog2, f,Writer=ascii.CommentedHeader)
        writefits(catalog2,mcats[1])

    catalog1 = catalog1._new_from_slice(imatch1)
    catalog2 = catalog2._new_from_slice(imatch2)
    
    if stack :
        outcat = table.hstack([catalog1,catalog2])
        outcat.add_column(distance)
        return outcat
    else :
        return (catalog1,catalog2)


def toRegionFile(catalog, filename, symbol = 'ellipse', subtag='',wcs=True):
#        
#   Dumps the catalog into a ds9 region file with symbols = "ellipse" or "point".
#   
    if subtag :
        subtag = '_'+subtag    
    f = open(filename,'w')
    writer = Writer(f)
    if wcs :
        coordtag='_WORLD'
        coordid='linear;'
    else :
        coordtag='_IMAGE'
        #coordid='image;'
        coordid='image;'
        
    fulltag = coordtag+subtag
    
    if symbol=="ellipse" :
        for i in xrange(0,len(catalog)) :
            writer.write_row([coordid, 'ellipse',
                              catalog['X'+fulltag][i],
                              catalog['Y'+fulltag][i],
                              catalog['A'+fulltag][i],
                              catalog['B'+fulltag][i],
                              catalog['THETA'+fulltag][i]])
    if symbol=="point" :
        for i in xrange(0,len(catalog)) :
            writer.write_row([coordid, 'x point ',
                              catalog['X'+fulltag][i],
                              catalog['Y'+fulltag][i]])
    f.close()

def writefits(table,file):
    fitstable=fits.HDUList([fits.PrimaryHDU(),fits.BinTableHDU.from_columns(table.as_array())])
    fitstable.writeto(file,clobber=True)

def readfits(catfile, imgext=None):
    fitscat=fits.open(catfile)
    data = []
    if imgext is None:
        ext = slice(1, None)
    else:
        ext = slice(2*imgext, 2*imgext + 1)
    for hdu in fitscat[ext] :
        if hdu.header['EXTNAME'] == 'LDAC_OBJECTS':
            data.append(hdu.data)
    return table.Table(np.hstack(data))

def read(catfile,format=None):
    if format=="fits":
        return readfits(catfile)
    elif format=="ext" :
        return readext(catfile)
    else : 
        return ascii.read(catfile,format=format)

def readext(catfile):
    reader = ascii.CommentedHeader()
    with open(catfile,'r') as f:
        header = f.readline()
        try :
            while header[0] != '#' and header != '':
                header = f.readline()
        except IndexError:
            logger.error('No header found in catalog file...')
            sys.exit()
            
        header = header[1:-1].split()
        ncols = len(header)

    reader.data.splitter.process_line = lambda x:process_line(x,ncols)
    return reader.read(catfile)
    
def process_line(x,ncols):
    splitx=x.split()
    nx = len(splitx)
    deltacol = ncols - nx
    if deltacol > 1:
        splitx += ['0']*deltacol
    
    return ' '.join(splitx)
    

def mag_zeropoint_ccd(zp_mag_ccd, list_ext_num, idx_ext):
    '''Returns zero point magnitude of ccd corresponding to index of extension
    of fits catalog
    - Parameters:
        + zp_mag_ccd: array of zeropoint magnitude of all CCD
        + list_ext_num: list of extension number already in order by KIDS CCD
        + idx_ext: the extension number: 1st, 2nd, ...
    '''
    mag_zeropoint = zp_mag_ccd[list_ext_num.index(idx_ext)]
    return mag_zeropoint


def correct_mag_ccd(catalog, zp_mag_ccd, list_ext_num, idx_ext):
    '''Correct the magnitude by the zero point magnitude of CCD
    - Default magnitude: mag = -2.5*log10(flux) + MAG_ZEROPOINT
    - Corrected magnitude: mag_corr = mag - MAG_ZEROPOINT + MAG_ZEROPOINT_CCD
    Parameters:
        - catalog: input fits catalog
        - zp_mag_ccd: zero point magnitude of all ccd
        - list_ext_num: list of extension numbers corresponds to the CCD number
        - idx_ext: index of extension of catalog to correct magnitude
    Returns:
        - array of corrected magnitude
    '''
    cat_data = catalog[2*idx_ext].data
    mag = cat_data['MAG_AUTO']

    MAG_ZEROPOINT_CCD = mag_zeropoint_ccd(zp_mag_ccd, list_ext_num, idx_ext)
    mag_corr = mag - MAG_ZEROPOINT + MAG_ZEROPOINT_CCD
    return mag_corr


class Writer :
    def __init__(self, f):
        self.f = f
        
    def close(self):
        self.f.close()
        
    def write_comment(self, comment):
        self.f.write("# "+comment+NEWLINE)
        
    def write_row(self, row):
        self.f.write(" ".join([str(e) for e in row])+NEWLINE)


        
