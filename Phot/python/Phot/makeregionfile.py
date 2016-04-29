"""
@file: python/Phot/makeregionfile.py
@author: user
@date: 03/24/16
"""

import argparse
import ElementsKernel.Logging as log
from astropy.io import ascii
from . import catalog,utils

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
    parser.add_argument('--wcs', action='store_true',
                      help='Region file in wcs coordinates',
                      default=False)
    parser.add_argument('--catalog', type=str, help='Ascii catalog file.')
    parser.add_argument('--format',  type=str, help='Format of input catalog')
    parser.add_argument('--symbol',  type=str, help='Symbol to draw')
    parser.add_argument('--xkey',    type=str, help='Keyword for x position')
    parser.add_argument('--ykey',    type=str, help='Keyword for y position')
    parser.add_argument('--akey',    type=str, help='Keyword for ellipse long axis')
    parser.add_argument('--bkey',    type=str, help='Keyword for ellipse short axis')
    parser.add_argument('--thetakey',type=str, help='Keyword for ellipse orientation')
    parser.add_argument('--subtag',  type=str, help='subtag for subcatalog designation (ex : column X_WORLD would become X_WORLD_subtag)')
    #

    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is 
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger('makeregionfile')

    logger.info('#')
    logger.info('# Entering makeregionfile mainMethod()')
    logger.info('#')

    # !! Getting the option from the example option in defineSpecificProgramOption 
    # !! e.g string_option = args.string_value

    #
    cat = ascii.read(args.catalog,format=args.format)
    regionfile = utils.rm_extension(args.catalog)+'.reg'
    if args.subtag is None :
        args.subtag = ''
        
    catalog.toRegionFile(cat, regionfile, symbol = args.symbol, subtag=args.subtag, wcs=args.wcs)
    #

    logger.info('#')
    logger.info('# Exiting makeregionfile mainMethod()')
    logger.info('#')
