"""
@file python/Phot/photometry.py
@date 03/16/16
@author user
"""
from __future__ import division
import ElementsKernel.Logging as log
logger = log.getLogger('image')
import subprocess, os, sys
from . import utils
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
import numpy as np
import numpy.ma as ma
from string import Template
from matplotlib import pylab as P
import json
from scipy.ndimage.measurements import center_of_mass


def parse_section_list(sectionlist):
    section = [sys.maxint, 0, sys.maxint, None]
    for sectionstring in sectionlist :
        s = parse_section_string(sectionstring)
        section = [min(s[0], section[0]), max(s[1],section[1]), min(s[2],section[2]), max(s[3],section[3])]
    return section

def parse_section_string(section):
    section=section[1:-1].replace(':',',').split(',')
    section=map(int, section)
    return section

def parse_section(section):
    if type(section) is str:
        return parse_section_string(section)
    elif type(section) is list:
        return parse_section_list(section)

def get_section(header, *sections, **kwargs):
    section = [sys.maxint, 0, sys.maxint, None]
    for s in sections :
        if 'i' in kwargs and kwargs['i'] != '' :
            s = parse_section(header[s][uppercase.index(kwargs['i'])])
        else:
            s = parse_section(header[s])
        section = [min(s[0], section[0]), max(s[1],section[1]), min(s[2],section[2]), max(s[3], section[3])]
    return section

def get_pixscale(wcsdict):
    return np.sqrt(abs(np.linalg.det([[wcsdict['CD1_1'], wcsdict['CD1_2']], [wcsdict['CD2_1'], wcsdict['CD2_2']]])))*3600.

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

def plotfits(imname, ext=1, show=False, **kwargs):
    hdu = fits.open(imname)[ext]
    wcs = get_wcs(hdu.header)
    fig = P.figure()
    fig.add_subplot(111, projection=wcs)
    P.imshow(hdu.data, origin='lower', norm=LogNorm(), cmap='gray_r',**kwargs)
    P.colorbar()
    P.xlabel('RA')
    P.ylabel('Dec')

    return wcs

def plotcat(cat, wcs, show=False, poskeys=['X_WORLD', 'Y_WORLD'], **kwargs):
    objects_pos = wcs.all_world2pix(cat[poskeys[0]], cat[poskeys[1]], 0)
    P.scatter(objects_pos[0], objects_pos[1], **kwargs)
    if show :
        P.show()
    return objects_pos

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def mag2adu(mag,zp,exptime=1.):
    aduflux = exptime*10.**((zp-mag)/2.5)
    return aduflux

def smag2pix(smag,pixscale,zp,exptime=1.):
    pixvalue = exptime*10.**((zp-smag)/2.5)*pixscale**2
    return pixvalue

def sigma_background(counts_bg, readout_noise, gain):
    return np.sqrt(1/gain*(counts_bg + readout_noise))

def read_instrument(instrument_file):
    with open(instrument_file) as f :
        instrument = json.load(f)
    instrument["CCDS_LAYOUT"] = {int(k):v for (k, v) in instrument["CCDS_LAYOUT"].items()}
    return instrument

def updatewcs(im,instrument):
    for i,ext in enumerate(im[1:],1):
        wcs=instrument['CCDS_LAYOUT'][i]["WCS"]
        for k,v in wcs.items():
            ext.header[k]=v

def build_circle_mask(radius, pos=(0,0)):
    intradius = int(np.ceil(radius))
    y,x = np.ogrid[-intradius: intradius+1, -intradius: intradius+1]
    boolmask = x**2+y**2 <= radius**2
    mask = np.where(boolmask == True)
    mask = (mask[0] - intradius + pos[1], mask[1] - intradius + pos[0])
    return mask


def build_square_mask(radius, pos=(0,0)):
    intradius=int(np.ceil(radius))
    mask = np.mgrid[-intradius: intradius+1, -intradius: intradius+1]
    mask = (mask[0] + pos[1], mask[1] + pos[0])
    return mask

def measure_psfs_at(im, positions, radius, pixscale=1., show=False):
    psfs = []
    for pos in positions:
        mask = build_square_mask(radius, pos=pos)
        p = measure_psf(im.data[mask], pixscale=pixscale, show=show)
        psfs.append(2.3548*np.mean([p.x_stddev.value, p.y_stddev.value]))
    return np.array(psfs)

def measure_psfs(vignets, pixscale=1, show=False, mask_value=-1e30):
    psfs=[]
    for vignet in vignets :
        p = measure_psf(vignet, pixscale=pixscale, mask_value=mask_value, show=show)
        psfs.append(2.3548*np.mean([p.x_stddev.value, p.y_stddev.value]))
    return np.array(psfs)

def measure_barycenter(vignet, pixscale=1.):
    x,y = center_of_mass(vignet)
    x -= (vignet.shape[0]-1)/2
    x *= pixscale
    y -= (vignet.shape[1]-1)/2
    y *= pixscale    
    return (x,y)
    
def measure_psf(vignet, pixscale=1., show=False, mask_value=None):
    y, x = np.mgrid[-vignet.shape[0]/2:vignet.shape[0]/2, -vignet.shape[1]/2:vignet.shape[1]/2]*pixscale
    if mask_value :
        vignet = ma.masked_values(vignet, mask_value).filled(0)
    # Fit the data using astropy.modeling
    p_init=models.Gaussian2D(amplitude=vignet.max(), x_mean=0., y_mean=0.,
        x_stddev=2*pixscale, y_stddev=2*pixscale, theta=0, cov_matrix=None)
    fit_p = fitting.LevMarLSQFitter()

    p = fit_p(p_init, x, y, vignet)
    barycenter=measure_barycenter(vignet, pixscale=pixscale)
    
    # Plot the data with the best-fit model
    P.figure(figsize=(8, 2.5))
    P.subplot(1, 3, 1)
    P.imshow(vignet, origin='lower', interpolation='nearest', vmin=vignet.min(), vmax=vignet.max())
    P.title("Data")
    P.subplot(1, 3, 2)
    P.imshow(p(x, y), origin='lower', interpolation='nearest', vmin=vignet.min(), vmax=vignet.max())
    P.scatter(vignet.shape[0]/2, vignet.shape[1]/2,marker="+")
    P.annotate("({:.3f},{:.3f})".format(*barycenter), (vignet.shape[0]/3, vignet.shape[1]/3))
    P.title("Model - psf = {:.2f}".format(2.3548*np.mean([p.x_stddev.value, p.y_stddev.value])))
    P.subplot(1, 3, 3)
    P.imshow(vignet - p(x, y), origin='lower', interpolation='nearest', vmin=-vignet.max()/10,vmax=vignet.max()/10)
    P.title("Residual")
    P.tight_layout()
    if show :
        P.show()
    
    return p
    
def simple_aper_phot(im,positions,radius):
    fluxes=[]
    for pos in positions :
        mask=build_circle_mask(radius, pos=pos)
        fluxes.append(sum(im.data[mask]))
    return np.array(fluxes)

def build_annulus_mask(in_radius, out_radius, pos=(0,0)):
    ''' Build annulus mask to extract the flux inside a ring shape from the
    center of the object
    Parameters:
        - in_radius: inside radius (smaller distance from the center)
        - out_radius: outside radius
        - pos: position of the object
    '''
    intradius=int(np.ceil(out_radius))
    y,x = np.ogrid[-intradius: intradius+1, -intradius: intradius+1]
    boolmask = (x**2+y**2 >= in_radius**2) & (x**2+y**2 <= out_radius**2)
    mask = np.where(boolmask==True)
    mask = (mask[0] - intradius + pos[1], mask[1] - intradius + pos[0])
    return mask

def annulus_photometry(im, positions, in_radius, out_radius):
    ''' Intergral flux of a ring shape from the center of an object inside
        an image
    Parameters:
        - im: image data (n*m pixel)
        - positions: center position of object
        - in_radius: inside radius (smaller distance from the center)
        - out_radius: outside radius
    '''
    fluxes = []
    for pos in positions:
        mask = build_annulus_mask(in_radius, out_radius, pos=pos)
        fluxes.append(sum(im.data[mask]))
    #total_bg_pixels = np.shape(mask[1])
    return np.array(fluxes)


def sex(imname, zeropoint=0, outputcat=None):
    with open(utils.getAuxPathFile("default.sex.template"),'r') as f :
        conftemp = Template(f.read())
    conf = conftemp.substitute({'zeropoint':str(zeropoint)})
    with open("default.sex", 'w') as f :
        f.write(conf)
    if outputcat is None :
        outputcat = utils.rm_extension(imname)+".cat"

    cmd = ["sex",imname,"-c", "default.sex"]
    cmd += ["-CATALOG_NAME", outputcat]
    cmd += ["-PARAMETERS_NAME", utils.getAuxPathFile("default.param")]
    cmd += ["-FILTER_NAME", utils.getAuxPathFile("gauss_5.0_9x9.conv")]
    cmd += ["-STARNNW_NAME", utils.getAuxPathFile("default.nnw")]

    p = subprocess.call(cmd)

    if p!=0 :
        sys.exit("SExtractor failed... Exiting.")

def generate_weight_map(instrument):
    wmap = fits.HDUList([fits.PrimaryHDU()])
    for ccd in instrument["CCDS_LAYOUT"].values():
        data = np.zeros((ccd["HEIGHT"], ccd["WIDTH"]))
        datasec = get_section(ccd, "DATASEC")
        data[datasec[2]-1:datasec[3]-1, datasec[0]-1:datasec[1]-1]=1
        wmap.append(fits.ImageHDU(data.astype('uint8')))
    return wmap

def get_zeropoints(imname, apply_exptime=False, zerokey="SIMMAGZP"):
    im = fits.open(imname)
    zeropoints=[]
    for ext in im[1:]:
        z=ext.header[zerokey]
        if apply_exptime :
            z += 2.5*np.log10(im[0].header['EXPTIME'])
        zeropoints.append(z)

    logger.info("Mag zeropoints : {}".format(zeropoints))
    return zeropoints

def flux2mag(flux, zeropoint, exptime):
    return -2.5*np.log10(flux*1./exptime) + zeropoint

def get_ccd_wcs(telescope_pointing, instrument):
    ccd = instrument['CCDS_LAYOUT'][1]
    ccd['WCS'][u'CRVAL1']  = telescope_pointing['ra']
    ccd['WCS'][u'CRVAL2'] = telescope_pointing['dec']
    wcs_ccd = WCS(ccd['WCS'])
    return wcs_ccd

def generate_fake_star_catalog(telescope_pointing, instrument, step, objtype, magnitude):
    ccd = instrument['CCDS_LAYOUT'][1]
    wcs_ccd = get_ccd_wcs(telescope_pointing, instrument)
    with open('my_star_pixel_catalog.txt', 'w') as f_pix:
        with open("my_star_catalog.txt", 'w') as f_world:
            f_pix.write("# Objtype X_IMAGE Y_IMAGE \n")
            f_world.write("# X_WORLD Y_WORLD MAG \n")
            for i in xrange(step, ccd['WIDTH'], step):
                for j in xrange(step, ccd['HEIGHT'], step):
                    f_pix.write("%d %d %d" %(objtype, i, j) + "\n")
                    ra, dec = wcs_ccd.all_pix2world(i, j, 1)
                    f_world.write("%d %f %f %.2f" %(objtype, ra, dec, magnitude) + "\n")

def weighted_center(im):
    x = 0
    y = 0
    weight = 0
    for i in xrange(np.shape(im)[0]):
        for j in xrange(np.shape(im)[1]):
            x += (i+1)*im[i,j]
            y += (j+1)*im[i,j]
            weight += im[i,j]
    x_center = x/weight
    y_center = y/weight
    return x_center, y_center

def weighted_std(im, x_center, y_center):
    x_std2 = 0
    y_std2 = 0
    weight = 0
    for i in xrange(np.shape(im)[0]):
        for j in xrange(np.shape(im)[1]):
            x_std2 += ((i+1)-x_center)**2
            y_std2 += ((j+1)-y_center)**2
            weight += im[i,j]
    x_std = x_std2/weight
    y_std = y_std2/weight
    return x_std, y_std2
    