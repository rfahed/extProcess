"""
@file tests/python/smoke_test.py
@date 12/16/15
@author user
"""
import py.test
import os.path
import glob
from extSim import extsim, utils
import shutil
from fixtures.PyDataSyncFixture import *

data_dir = os.path.join(os.environ["PHOTROOT"],"tests", "data")

def symlinks(datafiles,workspace):
    for dst,src in datafiles.items():
        os.symlink(os.path.join(data_dir,src), os.path.join(workspace,dst))

def simulate():
    datafiles={"catalog.txt":"myEXTrandomcat.txt", "target.json":"des_target.json", "instrument.json":"des_test_instrument.json", "flat.fits":"des_testflat.fits", "bias.fits":"des_testbias.fits", "stamps.fits":"EUC-TEST-STAMPS-2017-02-21T154652.049159.fits"}
    with utils.mock_workspace('test_smoke_stamps_ws_',del_tmp=False) as workspace:
       symlinks(datafiles,workspace)
       args = extsim.parseOptions(['--workspace',workspace],defaultconfig='smoke_test.conf')
       extsim.mainMethod(args)
    return args

class Teststamps(object):
    """
    @class Testsmoke
    @brief Unit Test class
    """

    def setup_class(self): 
        PyDataSyncFixture(os.path.join(data_dir,"config/sync.conf"), os.path.join(data_dir,"config/test_file_list.txt"))
        self.del_tmp = True
        self.silent = True
        self.args = simulate()
        

    def test_outfile(self):
        assert glob.glob(os.path.join(self.args.output_path, '*.fits'))

    def teardown_class(self):
        shutil.rmtree(self.args.workspace)