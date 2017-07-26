"""
@file python/Phot/utils.py
@date 03/16/16
@author user
"""

import re
import os
import shutil

def make_figures_dir(test_name=""):
    try :
        basedir = os.environ["EXTSIM_FIGURES_DIR"]
    except KeyError :
        basedir = os.environ["PWD"]
    
    figdir = os.path.join(basedir,"extSimFigures",test_name)
    if not os.path.exists(figdir):
        os.makedirs(figdir)
    else:
        shutil.rmtree(figdir)           
    os.makedirs(figdir)
    return figdir


def rm_extension(filename):
    matches = re.search('(.+)\..+',filename)
    if matches :
        return matches.group(1)
    else :
        return filename
    
def getAuxPathFile(file_name):
    """
    Look for the <auxdir> path in the <ELEMENTS_AUX_PATH> environment variable 
    where is located the <auxdir/file_name> file. It returns the filename with 
    the path or an empty string if not found.
    """
    found = False
    full_filename = ''
    aux_dir = os.environ.get('ELEMENTS_AUX_PATH')
    if not aux_dir is None:
        for elt in aux_dir.split(os.pathsep):
            # look for the first valid path
            full_filename = os.path.sep.join([elt, file_name])
            if os.path.exists(full_filename) and 'auxdir' in full_filename:
                found = True
                break

    if not found:
        logger.error("# Auxiliary file NOT FOUND  : <%s>" % full_filename)
        full_filename = ''

    return full_filename