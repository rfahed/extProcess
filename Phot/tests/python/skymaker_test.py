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
    with utils.mock_workspace('test_skymaker_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extsim.parseOptions(['--workspace',workspace, '--truth_source', 'FAKESTARS'],defaultconfig='skymaker_test.conf')
       extsim.mainMethod(args)
    return args

class Testskymaker(object):
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
        
        catalog.scattercols(mergedcat,'DELTAX','DELTAY',xlab='Delta X (pixels)',ylab='Delta Y (pixels)',show=True)
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("test_positions_1.png")
        
        p.figure()
        f, axarr = p.subplots(2, sharex=True)
        axarr[0].hist(mergedcat['DELTAX'],bins=np.linspace(-3,3,61))
        axarr[0].set_title('Delta X (pixels)')
        axarr[1].hist(mergedcat['DELTAY'],bins=np.linspace(-3,3,61))
        axarr[1].set_title('Delta Y (pixels)')
        
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("test_positions_2.png")
        
        assert True
        
    def teardown_class(self):
        """
        Removes workspace
        """
        #shutil.rmtree(self.args.workspace, ignore_errors=True)
        

        
                    