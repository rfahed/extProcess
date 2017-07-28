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


data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={"target.json":"kids_target_nobg.json", "instrument.json":"oneccd_test_instrument.json"}
    with extSim.utils.mock_workspace('test_skymaker_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extSim.extsim.parseOptions(['--workspace',workspace],defaultconfig='skymaker_test.conf')
       extSim.extsim.mainMethod(args)
    return args

class Testskymaker(object):
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

    def test_mags_sex(self):
        tolflux = 10.
        tolmag = 0.01
        image.sex(os.path.join(self.args.output_path, 'output.fits'), zeropoint=self.target["zeropoint"]+2.5*np.log10(self.target["exptime"]))
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        outputcat=catalog.readfits(os.path.join(self.args.output_path, 'output.cat'))
        out_flux = np.median(outputcat["FLUX_AUTO"])
        out_mag = np.median(outputcat["MAG_AUTO"])
        expected_flux = image.mag2adu(inputcat["MAG"][0],self.target["zeropoint"],exptime=self.target["exptime"])
        expected_mag = inputcat["MAG"][0]
        assert np.abs(out_mag-expected_mag) < tolmag

    def test_mags(self):
        tol = 10.
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        pixradius = 5.*self.target["psf"]/self.instrument["PIXEL_SCALE"]
        positions = zip(inputcat["X_IMAGE"],inputcat["Y_IMAGE"])
        fluxes = image.simple_aper_phot(self.im[1],positions,pixradius)
        mean_flux = np.median(fluxes)
        expected_flux = image.mag2adu(inputcat["MAG"][0],self.target["zeropoint"],exptime=self.target["exptime"])
        assert np.abs(mean_flux-expected_flux) < tol

    def test_positions(self):
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        inputcat=catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
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
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,101))
        axarr[1].set_title('Delta Y (pixels)')

        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig(os.path.join(self.figdir,"positions_2.png"))
        tol = 0.05
        assert (np.mean(mergedcat['DELTAX']) < tol) and (np.mean(mergedcat['DELTAY']) < tol)

    def teardown_class(self):
        """
        Removes workspace
        """
        #shutil.rmtree(self.args.workspace, ignore_errors=True)
