"""
@file tests/python/smoke_test.py
@date 12/16/15
@author user
"""
import py.test
import os.path
import glob
import json
import extSim.utils
import extSim.extsim
from Phot import image, catalog, utils
from astropy.io import fits
import shutil
from string import uppercase as amptags
from fixtures.PyDataSyncFixture import *
from matplotlib import pylab as p
import subprocess
import numpy as np
from scipy.optimize import curve_fit
PSFFILE = "psf_1.3_pixscale-0.05.fits"

data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={"target.json":"des_target_nobg.json", "instrument.json":"des_oneccd.json", "psf_1A.fits":PSFFILE, "psf_1B.fits":PSFFILE, "psf.fits":PSFFILE, "stamps.fits":"stamps_gauss_pixscale-0.05.fits"}
    with extSim.utils.mock_workspace('test_stamp_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extSim.extsim.parseOptions(['--workspace',workspace],defaultconfig='stamp_test.conf')
       extSim.extsim.mainMethod(args)
    return args

def get_stamp_size(stampfile):
    stamps=fits.open(stampfile)
    pixscale=image.get_pixscale(stamps[0].header)
    p = image.measure_psf(stamps[0].data, pixscale=pixscale)
    return 2.3548*np.mean([p.x_stddev.value,p.y_stddev.value])
    
def get_psf_size(psffile):
    psfgrid=fits.open(psffile)
    i=int(psfgrid[0].data.shape[0]/2)
    j=int(psfgrid[0].data.shape[1]/2)
    vignet=psfgrid[0].data[i,j]
    p = image.measure_psf(vignet, pixscale=psfgrid[0].header["PSF_SAMP"])
    return 2.3548*np.mean([p.x_stddev.value,p.y_stddev.value])

class Teststamp(object):
    """
    @class Testsmoke
    @brief Unit Test class
    """

    def setup_class(self):
        PyDataSyncFixture("../../config/sync.conf", "../../config/test_file_list.txt")
        self.del_tmp = True
        self.silent = True
        self.args = simulate()
        self.figdir = utils.make_figures_dir(__name__)
        self.instrument = extSim.utils.read_instrument(os.path.join(self.args.workspace,self.args.instrument))
        with open(os.path.join(self.args.workspace,self.args.target)) as f:
            self.target = json.load(f)
        try :
            self.im = fits.open(os.path.join(self.args.output_path, 'output.fits'))
        except :
            self.im = None
        self.inputpsf = get_psf_size(os.path.join(self.args.workspace,"psf.fits"))*self.instrument["PIXEL_SCALE"]
        self.inputstampsize = get_stamp_size(os.path.join(self.args.workspace,"stamps.fits"))
        self.expected_outstampsize = np.sqrt(self.inputpsf**2+self.inputstampsize**2)

    def test_outfile(self):
        """
        Check that an output image have been generated in the output directory.
        """
        assert glob.glob(os.path.join(self.args.output_path, 'output.fits'))

    def test_mags_sex(self):
        tolflux = 10.
        tolmag = 0.01
        image.sex(os.path.join(self.args.output_path, 'output.fits'), zeropoint=self.target["zeropoint"]+2.5*np.log10(self.target["exptime"]))
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1B.cat'))
        outputcat=catalog.readfits(os.path.join(self.args.output_path, 'output.cat'))
        out_flux = np.median(outputcat["FLUX_AUTO"])
        out_mag = np.median(outputcat["MAG_AUTO"])
        expected_flux = image.mag2adu(inputcat["MAG"][0],self.target["zeropoint"],exptime=self.target["exptime"])
        expected_mag = inputcat["MAG"][0]
        assert np.abs(out_mag-expected_mag) < tolmag

    def test_mags(self):
        tol = 10.
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1B.cat'))
        pixradius = 5.*self.inputpsf/self.instrument["PIXEL_SCALE"]
        positions = zip(inputcat["X_IMAGE"],inputcat["Y_IMAGE"])
        fluxes = image.simple_aper_phot(self.im[1],positions,pixradius)
        mean_flux = np.median(fluxes)
        expected_flux = image.mag2adu(inputcat["MAG"][0],self.target["zeropoint"],exptime=self.target["exptime"])
        assert np.abs(mean_flux-expected_flux) < tol

    def test_positions(self):
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        inputcat=catalog.read(os.path.join(self.args.tmp_path, 'ccd_1B.cat'))
        outputcat=catalog.readfits(os.path.join(self.args.output_path, 'output.cat'))
        ouputcat=outputcat[outputcat["FLAGS"]==0]
        mergedcat=catalog.mergecatpix(inputcat,outputcat, delta=3, poskeys1=['X_IMAGE','Y_IMAGE'], poskeys2=['X_IMAGE','Y_IMAGE'])
        catalog.add_diff_column(mergedcat,'X_IMAGE_2','X_IMAGE_1',outputfield='DELTAX')
        catalog.add_diff_column(mergedcat,'Y_IMAGE_2','Y_IMAGE_1',outputfield='DELTAY')

        #p.figure()
        #p.scatter(mergedcat['DELTAX'],mergedcat['DELTAY'])
        #p.xlabel('Delta X (pixels)')
        #p.ylabel('Delta Y (pixels)')

        catalog.scattercols(mergedcat,'DELTAX','DELTAY',xlab='Delta X (pixels)',ylab='Delta Y (pixels)',show=False)
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig(os.path.join(self.figdir,"positions_1.png"))

        p.figure()
        f, axarr = p.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,101))
        axarr[0].set_title('Delta X (pixels)')
        p.grid()
        p.legend()
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,101))
        axarr[1].set_title('Delta Y (pixels)')
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig(os.path.join(self.figdir,"positions_2.png"))
        tol = 0.15
        assert (np.mean(mergedcat['DELTAX']) < tol) and (np.mean(mergedcat['DELTAY']) < tol)
    
    def test_stamp_size(self):
        tol=0.07
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1B.cat'))
        pixradius = 2.*self.inputpsf/self.instrument["PIXEL_SCALE"]
        positions = zip(inputcat["X_IMAGE"],inputcat["Y_IMAGE"])
        sizes=image.measure_psfs_at(self.im[1],positions,pixradius,pixscale=self.instrument["PIXEL_SCALE"], show=False)
        assert np.abs(np.median(sizes) - self.expected_outstampsize) < tol
    
    def teardown_class(self):
        """
        Removes workspace
        """
        #shutil.rmtree(self.args.workspace, ignore_errors=True)
