"""
@file: python/Phot/sextractor.py
@author: user
@date: 03/16/16
"""

import argparse
import ElementsKernel.Logging as log
from . import utils
from . import image, catalog

def defineSpecificProgramOptions():
    """
    @brief Allows to define the (command line and configuration file) options
    specific to this program
    
    @details
        See the Elements documentation for more details.
    @return
        An  ArgumentParser.
    """

    parser = argparse.ArgumentParser()

    #
    # !!! Write your program options here !!!
    parser.add_argument('image', metavar='image', type=str, help='input fits image')
    parser.add_argument('--zerokey', type=str, help='zeropoint key in fits header')
    parser.add_argument('--zeropoint', type=str, help='zeropoint type (PER_S or INTEGRATED)')
    parser.add_argument('--outputcat', type=str, help='Output catalog file name')
    
    #parser.add_argument('-c', '--clobber', action='store_true',
    #                  help='force overwriting files',
    #                  default=False,
    #                  dest='clobber')
    #

    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is 
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger('sextractor')

    logger.info('#')
    logger.info('# Entering sextractor mainMethod()')
    logger.info('#')
    
    if args.zeropoint == "PER_S":
        zeropoints = image.get_zeropoints(args.image,apply_exptime=True,zerokey=args.zerokey)
    elif args.zeropoint == "INTEGRATED":
        zeropoints = image.get_zeropoints(args.image,apply_exptime=False,zerokey=args.zerokey)
    
    if args.outputcat == None :
        args.outputcat=utils.rm_extension(args.image)+".cat"
    
    image.sex(args.image,outputcat=args.outputcat)
    
    catalog.apply_zeropoints(args.outputcat,zeropoints)
    
    # !! Getting the option from the example option in defineSpecificProgramOption 
    # !! e.g string_option = args.string_value

    
    
    logger.info('#')
    logger.info('# Exiting sextractor mainMethod()')
    logger.info('#')
