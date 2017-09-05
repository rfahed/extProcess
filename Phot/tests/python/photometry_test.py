"""
@file tests/python/photometry_test.py
@date 12/16/15
@author user
"""
import py.test
import os.path
import glob, json
import extSim.utils
import extSim.extsim
from Phot import image, catalog, utils
from astropy.io import fits
from astropy.wcs import WCS
from astropy import table
import shutil
from string import uppercase as amptags
from fixtures.PyDataSyncFixture import *
from matplotlib import pylab as p
import subprocess
import numpy as np
from scipy.optimize import curve_fit


data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={"target.json":"kids_target.json", "instrument.json":"oneccd_test_instrument.json"}
    with extSim.utils.mock_workspace('test_photometry_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extSim.extsim.parseOptions(['--workspace',workspace],defaultconfig='skymaker_test.conf')
       extSim.extsim.mainMethod(args)
    return args

class TestPhotometry(object):
    """
    @class Testsmoke
    @brief Unit Test class
    """

    def setup_class(self):
        PyDataSyncFixture("../../config/sync.conf", "../../config/test_file_list.txt")
        self.del_tmp = True
        self.silent = True
        self.args = simulate()
        self.instrument = extSim.utils.read_instrument(os.path.join(self.args.workspace,self.args.instrument))
        self.figdir = utils.make_figures_dir(__name__)
        with open(os.path.join(self.args.workspace,self.args.target)) as f:
            self.target = json.load(f)
        try :
            self.im = fits.open(os.path.join(self.args.output_path, 'output.fits'))
        except :
            self.im = None

    def test_outfile(self):
        """
        Check that an output image have been generated in the output directory.
        """
        assert glob.glob(os.path.join(self.args.output_path, 'output.fits'))

    def test_counts(self):
        '''
        Check the image magnitude as expected magnitude
        '''
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        counts = self.im[1].data.ravel()
        expected_bg_counts = image.smag2pix(self.im[1].header['SIMMAGBG'],
        	self.instrument['PIXEL_SCALE'], self.im[1].header['SIMMAGZP'],
        	exptime=self.im[0].header['EXPTIME'])

        p.figure()
        p.hist(counts)
        p.title("Expected counts = {:0.2f}".format(expected_bg_counts))
        p.savefig(os.path.join(self.figdir,"photometry_counts.png"))


    def test_flux(self):
        tol = 150.
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        pixradius = 3*self.target["psf"]/self.instrument["PIXEL_SCALE"]
        positions = list(zip(inputcat["X_IMAGE"]-1, inputcat["Y_IMAGE"]-1))
        fluxes = image.simple_aper_phot(self.im[1], positions, pixradius)
        sky_background = image.annulus_photometry(self.im[1], positions,
        	pixradius+5, pixradius+8)

        total_bg_pixels = np.shape(image.build_annulus_mask(pixradius+5, pixradius+8, positions[0]))[1]
        total_source_pixels = np.shape(image.build_circle_mask(pixradius,
        	positions[0]))[1]

        estimated_fluxes = fluxes - sky_background*1./total_bg_pixels*total_source_pixels

        estimated_magnitude = image.flux2mag(estimated_fluxes,
        	self.im[1].header['SIMMAGZP'], self.target["exptime"])

        expected_flux = image.mag2adu(17.5, self.target["zeropoint"][0],
        	exptime=self.target["exptime"])

        p.figure()
        p.hist(fluxes, bins=50)
        p.title('Expected flux: {:0.2f}, mean flux: {:1.2f}'.format(expected_flux, np.mean(estimated_fluxes)))
        p.savefig(os.path.join(self.figdir,'Fluxes.png'))

        assert np.all(np.abs(fluxes-expected_flux) < tol)

    def teardown_class(self):
        """
        Removes workspace
        """
        #shutil.rmtree(self.args.workspace, ignore_errors=True)
