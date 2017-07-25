"""
@file tests/python/smoke_test.py
@date 12/16/15
@author user
"""
import py.test
import os.path
import glob
import extSim.utils
import extSim.extsim
from astropy.io import fits
import shutil
from string import uppercase as amptags
from fixtures.PyDataSyncFixture import *
import ElementsKernel.Logging as logging

log = logging.getLogger('image')

data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={"catalog.txt":"myEXTrandomcat-sersics.txt", "target.json":"des_target.json", "instrument.json":"des_test_instrument.json", "flat.fits":"des_testflat.fits", "bias.fits":"des_testbias.fits"}
    with extSim.utils.mock_workspace('test_smoke_sersics_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extSim.extsim.parseOptions(['--workspace',workspace],defaultconfig='smoke_test.conf')
       extSim.extsim.mainMethod(args)
    return args

class Testsersic(object):
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
        try :
            self.im = fits.open(os.path.join(self.args.output_path, 'output.fits'))
        except :
            self.im = None
        
    def test_outfile(self):
        """
        Check that an output image have been generated in the output directory.
        """
        assert glob.glob(os.path.join(self.args.output_path, 'output.fits'))
        
    def test_amp_positions(self):
        """
        Check that the amps sections header keywords have the expected values.
        """
        for iext,ccdinfos in self.instrument["CCDS_LAYOUT"].items():
            for key in [k for k in ccdinfos if k.endswith("SEC")]:
                for iamp in range(self.instrument["AMPS_X_CCD"]):
                    amptag = amptags[iamp]
                    log.info("extension : {} amp : {} key : {}".format(iext,amptag,key))
                    assert self.im[iext].header[key+amptag]==ccdinfos[key][iamp]
    
    def test_ccd_order(self):
        """
        Check that the fits extensions are ordered as expected.
        """
        for iext,ccdinfos in self.instrument["CCDS_LAYOUT"].items():
            assert self.im[iext].header["CCDNUM"]==iext and  self.im[iext].header["EXTNAME"]==ccdinfos["EXTNAME"]
            
    def test_pixel_headkeys(self):
        """
        Check that the pixel related headkeys are coherent with the input values.
        """
        headkeys={"GAIN":"GAIN","SATURAT":"SATUR_LEVEL","RDNOISE":"READOUT_NOISE"}
        for iext,ccdinfos in self.instrument["CCDS_LAYOUT"].items():
            for hkey,pkey in headkeys.items():
                for iamp in range(self.instrument["AMPS_X_CCD"]):
                    amptag = amptags[iamp]
                    assert self.im[iext].header[hkey+amptag]==ccdinfos[pkey][iamp]
                    
    def test_extend_keyword(self):
        """
        Check that the pixel related headkeys are coherent with the input values.
        """
        assert self.im[0].header["EXTEND"]==True
    
    def test_obstype_keyword(self):
        """
        Check that OBSTYPE keyword complies with DES convention
        """
        assert self.im[0].header["OBSTYPE"]=="object"
        
    def teardown_class(self):
        """
        Removes workspace
        """
        shutil.rmtree(self.args.workspace, ignore_errors=True)
        
                    
