"""
@file tests/python/star_test.py
@date 09/06/17
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
import matplotlib.pyplot as plt
import subprocess
import numpy as np
from scipy.optimize import curve_fit
from astropy.wcs import WCS
from astropy import table

data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={ "catalog.txt":"my_star_catalog.txt",
                "target.json":"kids_target.json",
                "instrument.json":"oneccd_test_instrument.json"}
    with extSim.utils.mock_workspace('test_star_ws_',del_tmp=False) as workspace:
        symlinks(datafiles,workspace)
        args = extSim.extsim.parseOptions(['--workspace',workspace], defaultconfig='star_test.conf')
        extSim.extsim.mainMethod(args)
    return args

def get_psf_size(psffile):
    psfgrid=fits.open(psffile)
    i=int(psfgrid[0].data.shape[0]/2)
    j=int(psfgrid[0].data.shape[1]/2)
    vignet=psfgrid[0].data[i,j]
    # 
    p = image.measure_psf(vignet, pixscale=psfgrid[0].header["PSF_SAMP"], show=True)
    return 2.3548*np.mean([p.x_stddev.value,p.y_stddev.value])

class Testpsf(object):
    """
    @class Teststar
    @brief Unit Test class
    """
    def setup_class(self):
        PyDataSyncFixture("../../config/sync.conf", "../../config/test_file_list.txt")
        self.del_tmp = True
        self.silent = True
        self.args = simulate()
        self.figdir = utils.make_figures_dir(__name__)
        self.instrument = extSim.utils.read_instrument(
            os.path.join(self.args.workspace,self.args.instrument))
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

    def test_wcs_positions(self):
        inputcat=catalog.read(os.path.join(self.args.workspace, self.args.catalog))
        pixcat=catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        
        outputcat=catalog.pixel_to_world(pixcat, WCS(self.im[1].header),
            pos_keys=['X_IMAGE', 'Y_IMAGE'])
    
        mergedcat=catalog.mergecat(inputcat,table.Table(outputcat), delta=2e-4,
            poskeys1=['X_WORLD','Y_WORLD'], poskeys2=['X_WORLD','Y_WORLD'])
        
        catalog.add_diff_column(mergedcat, 'X_WORLD_2', 'X_WORLD_1', outputfield='DELTAX')
        catalog.add_diff_column(mergedcat, 'Y_WORLD_2', 'Y_WORLD_1', outputfield='DELTAY')
        mergedcat['DELTAX']*=3600.
        mergedcat['DELTAY']*=3600.
        
        plt.figure()
        plt.scatter(mergedcat['DELTAX'], mergedcat['DELTAY'])
        plt.xlabel('Delta X (arcsec)')
        plt.ylabel('Delta Y (arcsec)')      
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(self.figdir,"wcs_positions_1.png"))
        
        plt.figure()
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,61))
        axarr[0].set_title('Delta X (arcsec)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,61))
        axarr[1].set_title('Delta Y (arcsec)')
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(self.figdir,"wcs_positions_2.png"))
        tol = 0.0001
        assert np.all(np.abs(mergedcat['DELTAX']) < tol) and np.all(np.abs(mergedcat['DELTAY']) < tol)

    def test_positions_sextractor(self):
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        inputcat=catalog.read(os.path.join(self.args.workspace,self.args.catalog))

        #output catalog is extracted by SExtractor program
        outputcat=catalog.readfits(os.path.join(self.args.output_path, 'output.cat'))
        ouputcat=outputcat[(50 < outputcat["X_IMAGE"]) & (outputcat["X_IMAGE"] < self.im[1].header["NAXIS1"] - 50) ]
        ouputcat=outputcat[(50 < outputcat["Y_IMAGE"]) & (outputcat["Y_IMAGE"] < self.im[1].header["NAXIS2"] - 50) ]

        mergedcat=catalog.mergecat(inputcat, outputcat, delta=2e-4,
            poskeys1=['X_WORLD','Y_WORLD'], poskeys2=['X_WORLD','Y_WORLD'])
        catalog.add_diff_column(mergedcat, 'X_WORLD_2', 'X_WORLD_1', outputfield='DELTAX')
        catalog.add_diff_column(mergedcat, 'Y_WORLD_2', 'Y_WORLD_1', outputfield='DELTAY')
        mergedcat['DELTAX']*=3600.
        mergedcat['DELTAY']*=3600.
        
        plt.figure()
        plt.scatter(mergedcat['DELTAX'], mergedcat['DELTAY'])
        plt.xlabel('Delta X (arcsec)')
        plt.ylabel('Delta Y (arcsec)')      
        plt.grid()
        plt.tight_layout()
        plt.savefig(os.path.join(self.figdir,"positions_1.png"))
        
        plt.figure()
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,61))
        axarr[0].set_title('Delta X (arcsec)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,61))
        axarr[1].set_title('Delta Y (arcsec)')
        
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.figdir,"positions_2.png"))
        tol = 0.5
        assert (np.mean(mergedcat['DELTAX']) < tol) and (np.mean(mergedcat['DELTAY']) < tol)
  
    def test_pix_positions(self):
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        inputcat=catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        outputcat=catalog.readfits(os.path.join(self.args.output_path, 'output.cat'))
        ouputcat=outputcat[outputcat["FLAGS"]==0]
        mergedcat=catalog.mergecatpix(inputcat,outputcat, delta=3,
            poskeys1=['X_IMAGE','Y_IMAGE'], poskeys2=['X_IMAGE','Y_IMAGE'])

        catalog.add_diff_column(mergedcat, 'X_IMAGE_2', 'X_IMAGE_1', outputfield='DELTAX')
        catalog.add_diff_column(mergedcat, 'Y_IMAGE_2', 'Y_IMAGE_1', outputfield='DELTAY')
        
        catalog.scattercols(mergedcat,'DELTAX','DELTAY',xlab='Delta X (pixels)',ylab='Delta Y (pixels)',show=False)
        
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.figdir,"wcs_test_pix_positions_1.png"))
        
        plt.figure()
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,101))
        axarr[0].set_title('Delta X (pixels)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,101))
        axarr[1].set_title('Delta Y (pixels)')
        
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.figdir,"wcs_test_pix_positions_2.png"))
        tol = 0.15
        assert (np.mean(mergedcat['DELTAX']) < tol) and (np.mean(mergedcat['DELTAY']) < tol)

    def test_counts(self):
        '''
        Check the image magnitude as expected magnitude
        '''
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        counts = self.im[1].data.ravel()
        expected_bg_counts = image.smag2pix(self.im[1].header['SIMMAGBG'],
            self.instrument['PIXEL_SCALE'], self.im[1].header['SIMMAGZP'],
            exptime=self.im[0].header['EXPTIME'])

        plt.figure()
        plt.hist(counts)
        plt.title("Expected counts = {:0.2f}".format(expected_bg_counts))
        plt.savefig(os.path.join(self.figdir,"photometry_counts.png"))

    def test_flux(self):
        tol = 150.
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        pixradius = 3*self.target["psf"]/self.instrument["PIXEL_SCALE"]
        positions = list(zip(inputcat["X_IMAGE"].astype(int), inputcat["Y_IMAGE"].astype(int)))
        fluxes = image.simple_aper_phot(self.im[1], positions, pixradius)
        sky_background = image.annulus_photometry(self.im[1], positions,
            pixradius+5, pixradius+8)

        total_bg_pixels = np.shape(image.build_annulus_mask(pixradius+5, pixradius+8, positions[0]))[1]
        total_source_pixels = np.shape(image.build_circle_mask(pixradius,
            positions[0]))[1]

        estimated_fluxes = fluxes - sky_background*1./total_bg_pixels*total_source_pixels

        estimated_magnitude = image.flux2mag(estimated_fluxes,
            self.im[1].header['SIMMAGZP'], self.target["exptime"])

        expected_flux = image.mag2adu(17, self.target["zeropoint"][0],
            exptime=self.target["exptime"])

        plt.figure()
        plt.hist(fluxes, bins=50)
        plt.title('Expected flux: {:0.2f}, mean flux: {:1.2f}'.format(expected_flux, np.mean(estimated_fluxes)))
        plt.savefig(os.path.join(self.figdir,'Fluxes.png'))

        assert np.all(np.abs(fluxes-expected_flux) < tol)

    def test_psf_size(self):
        tol=0.05
        inputcat = catalog.read(os.path.join(self.args.tmp_path, 'ccd_1B.cat'))
        pixradius = 2.*self.inputpsf/self.instrument["PIXEL_SCALE"]
        positions = zip(inputcat["X_IMAGE"],inputcat["Y_IMAGE"])
        psfs=image.measure_psfs_at(self.im[1],positions,pixradius,pixscale=self.instrument["PIXEL_SCALE"],show=False)
        assert np.abs(np.median(psfs) - self.inputpsf) < tol

        
    def teardown_class(self):
        """
        Removes workspace
        """
        #shutil.rmtree(self.args.workspace, ignore_errors=True)
