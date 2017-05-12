"""
@file tests/python/smoke_test.py
@date 12/16/15
@author user
"""
import py.test
import os.path
import glob
from extSim import extsim, utils
from Phot import image, catalog 
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
    datafiles={"target.json":"kids_target.json", "instrument.json":"oneccd_test_instrument.json", "catalog.txt":"wcstest_catalog.txt"}
    with utils.mock_workspace('test_wcs_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extsim.parseOptions(['--workspace',workspace],defaultconfig='wcs_test.conf')
       extsim.mainMethod(args)
    return args

class Testwcs(object):
    """
    @class Testsmoke
    @brief Unit Test class
    """
    
    def setup_class(self): 
        PyDataSyncFixture(os.path.join(data_dir,"config/sync.conf"), os.path.join(data_dir,"config/test_file_list.txt"))
        self.del_tmp = True
        self.silent = True
        self.args = simulate()
        self.instrument = utils.read_instrument(os.path.join(self.args.workspace,self.args.instrument))
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
        inputcat=catalog.read(os.path.join(self.args.workspace,self.args.catalog))
        pixcat=catalog.read(os.path.join(self.args.tmp_path, 'ccd_1.cat'))
        
        outputcat=catalog.pixel_to_world(pixcat, WCS(self.im[1].header,), pos_keys=['X_IMAGE', 'Y_IMAGE'])
        
        mergedcat=catalog.mergecat(inputcat,table.Table(outputcat), delta=2e-4, poskeys1=['X_WORLD','Y_WORLD'], poskeys2=['X_WORLD','Y_WORLD'])
        
        catalog.add_diff_column(mergedcat,'X_WORLD_2','X_WORLD_1',outputfield='DELTAX')
        catalog.add_diff_column(mergedcat,'Y_WORLD_2','Y_WORLD_1',outputfield='DELTAY')
        mergedcat['DELTAX']*=3600.
        mergedcat['DELTAY']*=3600.
        #p.figure()
        #p.scatter(mergedcat['DELTAX'],mergedcat['DELTAY'])
        #p.xlabel('Delta X (pixels)')
        #p.ylabel('Delta Y (pixels)')
        
        catalog.scattercols(mergedcat,'DELTAX','DELTAY',xlab='Delta X (arcsec)',ylab='Delta Y (arcsec)',show=False)
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("wcs_test_wcs_positions_1.png")
        
        p.figure()
        f, axarr = p.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,61))
        axarr[0].set_title('Delta X (arcsec)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,61))
        axarr[1].set_title('Delta Y (arcsec)')
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("wcs_test_wcs_positions_2.png")
        tol = 0.0001
        assert np.all(np.abs(mergedcat['DELTAX']) < tol) and np.all(np.abs(mergedcat['DELTAY']) < tol)
        

           
    def test_positions(self):
        image.sex(os.path.join(self.args.output_path, 'output.fits'))
        inputcat=catalog.read(os.path.join(self.args.workspace,self.args.catalog))
        outputcat=catalog.readfits(os.path.join(self.args.output_path, 'output.cat'))
        ouputcat=outputcat[(50 < outputcat["X_IMAGE"]) & (outputcat["X_IMAGE"] < self.im[1].header["NAXIS1"] - 50) ]
        ouputcat=outputcat[(50 < outputcat["Y_IMAGE"]) & (outputcat["Y_IMAGE"] < self.im[1].header["NAXIS2"] - 50) ]
        mergedcat=catalog.mergecat(inputcat,outputcat, delta=2e-4, poskeys1=['X_WORLD','Y_WORLD'], poskeys2=['X_WORLD','Y_WORLD'])
        catalog.add_diff_column(mergedcat,'X_WORLD_2','X_WORLD_1',outputfield='DELTAX')
        catalog.add_diff_column(mergedcat,'Y_WORLD_2','Y_WORLD_1',outputfield='DELTAY')
        mergedcat['DELTAX']*=3600.
        mergedcat['DELTAY']*=3600.
        #p.figure()
        #p.scatter(mergedcat['DELTAX'],mergedcat['DELTAY'])
        #p.xlabel('Delta X (pixels)')
        #p.ylabel('Delta Y (pixels)')
        
        catalog.scattercols(mergedcat,'DELTAX','DELTAY',xlab='Delta X (arcsec)',ylab='Delta Y (arcsec)',show=False)
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("wcs_test_positions_1.png")
        
        p.figure()
        f, axarr = p.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,61))
        axarr[0].set_title('Delta X (arcsec)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,61))
        axarr[1].set_title('Delta Y (arcsec)')
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("wcs_test_positions_2.png")
        tol = 0.5
        assert np.all(np.abs(mergedcat['DELTAX']) < tol) and np.all(np.abs(mergedcat['DELTAY']) < tol)
  
    def test_pix_positions(self):
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
        p.savefig("wcs_test_pix_positions_1.png")
        
        p.figure()
        f, axarr = p.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-0.5,0.5,101))
        axarr[0].set_title('Delta X (pixels)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-0.5,0.5,101))
        axarr[1].set_title('Delta Y (pixels)')
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("wcs_test_pix_positions_2.png")
        tol = 0.15
        assert np.all(np.abs(mergedcat['DELTAX']) < tol) and np.all(np.abs(mergedcat['DELTAY']) < tol)
  
        
    def teardown_class(self):
        """
        Removes workspace
        """
        #shutil.rmtree(self.args.workspace, ignore_errors=True)
        

        
                    