"""
@file: python/Phot/mergecatalogs.py
@author: user
@date: 03/17/16
"""

from __future__ import division
import argparse
import ElementsKernel.Logging as log
from . import catalog,utils
from astropy.io import ascii

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
    parser.add_argument('catalogs', nargs='+', type=str,  help='Catalogs to merge')
    parser.add_argument('--outputcat', type=str, help='Output catalog')
    parser.add_argument('--filters', nargs='+', type=str, help='Filter names (used to tag the output catalogs)')
    parser.add_argument('--tol', type=float, help='Maximum distance for source association in degrees.')
    parser.add_argument('--format', nargs='+', type=str, help='Format of header')
    
    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is 
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger('mergecatalogs')

    logger.info('#')
    logger.info('# Entering mergecatalogs mainMethod()')
    logger.info('#')

    # !! Getting the option from the example option in defineSpecificProgramOption 
    # !! e.g string_option = args.string_value

    catalogs = []
    for i,cat in enumerate(args.catalogs) :
	if len(args.format) > 1:
		ft = args.format[i] 
        catalogs.append(catalog.read(cat,format=ft))
    
    mergedcat = catalog.mergecats(catalogs,delta=args.tol,filters=args.filters)
    with open(args.outputcat, 'w') as f :
        ascii.write(mergedcat, f,Writer=ascii.CommentedHeader)
        
    logger.info('# Percentage of cat 1 matched : {} %'.format(100*(len(mergedcat)/len(catalogs[0]))))
    logger.info('# Percentage of cat 2 matched : {} %'.format(100*(len(mergedcat)/len(catalogs[1]))))
    logger.info('# Size ratio cat 2 / cat 1 : {} %'.format(100*(len(catalogs[1])/len(catalogs[0]))))



    logger.info('#')
    logger.info('# Exiting mergecatalogs mainMethod()')
    logger.info('#')
