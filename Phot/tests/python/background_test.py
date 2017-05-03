"""
@file tests/python/smoke_test.py
@date 12/16/15
@author user
"""
import py.test
import os.path
import glob
from extSim import extsim, utils
from Phot import image 
from astropy.io import fits
import shutil
from string import uppercase as amptags
from fixtures.PyDataSyncFixture import *
from matplotlib import pylab as p

import numpy as np
from scipy.optimize import curve_fit


data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={"catalog.txt":"myEXTemptycat.txt", "target.json":"kids_target.json", "instrument.json":"oneccd_test_instrument.json"}
    with utils.mock_workspace('test_background_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extsim.parseOptions(['--workspace',workspace, '--flat_name', '', '--bias_name', ''],defaultconfig='smoke_test.conf')
       extsim.mainMethod(args)
    return args

class Testbackground(object):
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
           
    def test_fwhm(self):
        """
        Check that background fwhm is consistent with poisson + readout noise
        """
        nbins=60
        minb=0
        bins=np.linspace(minb,minb+nbins-1,nbins)
        n_sim, bins, patches = p.hist(self.im[1].data.ravel(), bins=bins, alpha=0.5,label="sim image")

        p.xlabel('pixel value (ADU)',size=20)
        p.ylabel('counts',size=20)
  
        bins = (bins[:-1] + bins[1:])/2

        i_s=np.where(n_sim==max(n_sim))
        i_s=i_s[0][0]
        p0_s=[5000,bins[i_s],10]

        p_s, var_matrix = curve_fit(image.gauss, bins, n_sim, p0=p0_s)

        gauss_s = image.gauss(bins,*p_s)
        p.plot(bins, gauss_s, '--',label='sigma = {:.1f}'.format(p_s[2]),linewidth=2)
        
        expected_bg_adu=image.smag2pix(self.im[1].header['SIMMAGBG'],self.instrument['PIXEL_SCALE'],self.im[1].header['SIMMAGZP'],exptime=self.im[0].header['EXPTIME'])
        expected_fwhm=image.sigma_background(expected_bg_adu, self.im[1].header['RDNOISE'], self.im[1].header['GAIN'])
        tol=1.
        
        p.title('expected fwhm = {:0.2f} ADUs'.format(expected_fwhm))
        p.grid()
        p.legend()
        p.tight_layout()
        p.savefig("test_fwhm.png")

        assert abs(p_s[2] - expected_fwhm) < tol
        
    def teardown_class(self):
        """
        Removes workspace
        """
        shutil.rmtree(self.args.workspace, ignore_errors=True)
        

        
                    