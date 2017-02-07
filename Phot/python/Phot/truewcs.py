#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#


"""
File: python/Phot/truewcs.py

Created on: 02/07/17
Author: centos
"""

from __future__ import division, print_function
from future_builtins import *

import argparse
import ElementsKernel.Logging as log
from . import image
from astropy.io import fits
import json

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
    parser.add_argument('images', nargs='+', type=str,  help='Images to update')
    parser.add_argument('--instrument', type=str, help='Instrument file to use for WCS information')
    #
    # !!! Write your program options here !!!
    # e.g. parser.add_argument('--string-value', type=str, help='A string option')
    #

    return parser


def mainMethod(args):
    """
    @brief The "main" method.
    @details
        This method is the entry point to the program. In this sense, it is
        similar to a main (and it is why it is called mainMethod()).
    """

    logger = log.getLogger('truewcs')

    logger.info('#')
    logger.info('# Entering truewcs mainMethod()')
    logger.info('#')
    instrument = image.read_instrument(args.instrument)
    for imfile in args.images :
        im=fits.open(imfile)
        image.updatewcs(im,instrument)
        im.writeto('truewcs_'+imfile)
        

    logger.info('#')
    logger.info('# Exiting truewcs mainMethod()')
    logger.info('#')
