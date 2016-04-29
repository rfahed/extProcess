"""
@file python/Phot/catalog.py
@date 03/18/16
@author user
"""

from astropy.coordinates import ICRS, SkyCoord
from astropy import table
from astropy import units as u
import numpy as np
import matplotlib.pyplot as p

NEWLINE = '\n'

def mergecats(catalogs=None,delta=1e-4,filters=None,poskeys=['X_WORLD','Y_WORLD'],stack=True):
    ncat = len(catalogs)
    if filters and stack :
        _poskeys = []
        for i in range(ncat) :
            catalogs[i] = tag_catalog(catalogs[i], filters[i])
            _poskeys.append([x+'_'+filters[i] for x in poskeys])
      
    mergedcat = catalogs[0]
    for i in range(1,ncat):
        if filters :
            mergedcat=mergecat(mergedcat,catalogs[i],poskeys1=_poskeys[0],poskeys2=_poskeys[i],delta=delta,stack=stack)
        else :
            mergedcat=mergecat(mergedcat,catalogs[i], poskeys1=poskeys, poskeys2=poskeys, delta=delta,stack=stack)
            
    return mergedcat

def add_diff_column(catalog,field1,field2,outputfield=None):
    if outputfield is None:
        outputfield = "diff_{}_{}".format(field1,field2)
    diff = table.Column(name=outputfield,data=catalog[field1]-catalog[field2])
    catalog.add_column(diff)
    
    
def plotcols(catalog,field1,field2,show=False):    
    p.plot(catalog[field1],catalog[field2])
    p.grid()
    if show:
        p.show()
        
def scattercols(catalog,field1,field2,show=False):    
    p.scatter(catalog[field1],catalog[field2],marker='+')
    p.grid()
    if show:
        p.show()
    
def histogramcol(catalog,field,nbins,show=False,**kwargs):
    n, bins, patches = p.hist(np.array(catalog[field]),50,**kwargs)
    if show:
        p.show()
    
def tag_catalog(catalog, tag):
    for cname in catalog.colnames :
        catalog.rename_column(cname,cname+'_'+tag)
    return catalog

def matchcats(p1,p2,delta=1e-4):
    index,dist2d,dist3d = p1.match_to_catalog_sky(p2)
    imatch = np.where(dist2d < delta*u.degree)
    return (imatch, index[imatch])

def mergecat(catalog1,catalog2, delta=1e-4, poskeys1=['X_WORLD','Y_WORLD'], poskeys2=['X_WORLD','Y_WORLD'], stack=True):    
    p1=SkyCoord(catalog1[poskeys1[0]],catalog1[poskeys1[1]],frame='icrs',unit="deg")
    p2=SkyCoord(catalog2[poskeys2[0]],catalog2[poskeys2[1]],frame='icrs',unit="deg")
    
    imatch1, imatch2 = matchcats(p1,p2,delta=delta)
    
    catalog1 = catalog1._new_from_slice(imatch1)
    catalog2 = catalog2._new_from_slice(imatch2)
    
    if stack :
        return table.hstack([catalog1,catalog2])
    else :
        return (catalog1,catalog2)


def toRegionFile(catalog, filename, symbol = 'ellipse', subtag='',wcs=False):
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
    
class Writer :
    def __init__(self, f):
        self.f = f
        
    def close(self):
        self.f.close()
        
    def write_comment(self, comment):
        self.f.write("# "+comment+NEWLINE)
        
    def write_row(self, row):
        self.f.write(" ".join([str(e) for e in row])+NEWLINE)


        
