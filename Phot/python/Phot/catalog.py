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
    

def correct_mag_ccd(catalog, idx_ext, t_exp):
    '''Correct the magnitude by the zero point magnitude of CCD
    - Default magnitude: mag = -2.5*log10(flux) + MAG_ZEROPOINT
    - Corrected magnitude: mag_corr = mag - MAG_ZEROPOINT + MAG_ZEROPOINT_CCD
    Parameters:
        - catalog: input fits catalog
        - zp_mag_ccd: zero point magnitude of all ccd
        - list_ext_num: list of extension numbers corresponds to the CCD number
        - idx_ext: index of extension of catalog to correct magnitude
        - t_exp: exposure time of filter (telescope)
    Returns:
        - array of corrected magnitude
    '''
    cat_header = catalog[2*idx_ext-1].data
    cat_data = catalog[2*idx_ext].data
    mag = cat_data['MAG_AUTO']
    ## Get value of SIMMAGZP in header of output catalog
    for i in xrange(55, 65):
        if fits.Card.fromstring(cat_header[0][0][i])[0] == 'SIMMAGZP':
            MAG_ZEROPOINT_CCD = fits.Card.fromstring(cat_header[0][0][i])[1] 
    mag_corr = mag + 2.5*np.log10(t_exp) + MAG_ZEROPOINT_CCD
    return mag_corr

def cn_PnPoly(P, V):
    cn = 0    # the crossing number counter

    # repeat the first vertex at end
    V = tuple(V[:])+(V[0],)

    # loop through all edges of the polygon
    for i in range(len(V)-1):   # edge from V[i] to V[i+1]
        if ((V[i][1] <= P[1] and V[i+1][1] > P[1])   # an upward crossing
            or (V[i][1] > P[1] and V[i+1][1] <= P[1])):  # a downward crossing
            # compute the actual edge-ray intersect x-coordinate
            vt = (P[1] - V[i][1]) / float(V[i+1][1] - V[i][1])
            if P[0] < V[i][0] + vt * (V[i+1][0] - V[i][0]): # P[0] < intersect
                cn += 1  # a valid crossing of y=P[1] right of P[0]

    return cn % 2   # 0 if even (out), and 1 if odd (in)

def inside_footprint(catalog, footprint, pos_keys=['X_WORLD', 'Y_WORLD']):
    x = catalog[pos_keys[0]]
    y = catalog[pos_keys[1]]
    in_out = [cn_PnPoly([x[i], y[i]], footprint) for i in range(len(catalog))]
    return catalog[np.asarray(in_out, dtype=bool)]


def object_inside_footprint(outside_cat, inside_cat,
    pos_out=['X_WORLD', 'Y_WORLD'], pos_in=['X_WORLD', 'Y_WORLD']):

    x_min_out = min(inside_cat[pos_in[0]])
    x_max_out = max(inside_cat[pos_in[0]])
    y_min_out = min(inside_cat[pos_in[1]])
    y_max_out = max(inside_cat[pos_in[1]])
    select1 = (outside_cat[pos_out[0]] <= x_max_out) & (outside_cat[pos_out[0]] >= x_min_out)
    select2 = (outside_cat[pos_out[1]] <= y_max_out) & (outside_cat[pos_out[1]] >= y_min_out)
    return outside_cat[[select1&select2]]

def modify_position(cat, wcs, type_pos='PIXEL'):
    x, y = wcs.all_world2pix(cat['X_WORLD'], cat['Y_WORLD'], 0)
    if type_pos == 'WORLD':
        x, y = wcs.all_pix2world(cat['X_IMAGE'], cat['Y_IMAGE'], 0)
    return x, y


def match_catalogs(cat1, cat2, sep=0.0001, poskeys1=['X_WORLD','Y_WORLD'],
    poskeys2=['X_WORLD','Y_WORLD']):
    '''Find objects appear in both catalogs
    Parameters
        - cat1, cat2: two catalogs to match
        - sep: The on-sky separation to search within. unit: degree
        - poskeys1, poskeys2: position keywords of two catalogs.
    Return index of match objects in both catalogs
    '''
    c1 = SkyCoord(ra=cat1[poskeys1[0]], dec=cat1[poskeys1[1]], unit="deg")
    c2 = SkyCoord(ra=cat2[poskeys2[0]], dec=cat2[poskeys2[1]], unit="deg")
    idxc1, idxc2, d2d, d3d = c2.search_around_sky(c1, sep*u.deg)
    return idxc1, idxc2, d2d


class Writer :
    def __init__(self, f):
        self.f = f
        
    def close(self):
        self.f.close()
        
    def write_comment(self, comment):
        self.f.write("# "+comment+NEWLINE)
        
    def write_row(self, row):
        self.f.write(" ".join([str(e) for e in row])+NEWLINE)


        
