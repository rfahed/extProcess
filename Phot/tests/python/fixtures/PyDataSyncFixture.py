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
File: python/CT_TDM_Module/PyDataSyncFixture.py

Created on: 01/31/17
Author: user
"""

from __future__ import division, print_function
from future_builtins import *

import argparse
import os
import subprocess
import ElementsKernel.Logging as logging

log = logging.getLogger('PyDataSyncFixture')

class DataHost:
    """The test data hosting solution.
    """
    
    IRODS = "iRODS"
    WEBDAV = "WebDAV"
    DSS = "DSS"
    UNKNOWN = "unknown"
    
    __all__ = [IRODS, WEBDAV, DSS, UNKNOWN]
    
    @staticmethod
    def is_valid(p):
        return p in __all__
    
    #TODO make iterable

class PyDataSyncFixture(object):
    """Retrieve test data from a data repository prior to the unit test run.
    """
    
    ##########
    ## INIT ##
    ##########
    
    def __init__(self, connection=None, dependencies=None, host=None):
        """Create a fixture with given configuration.
        Three calls are possible:
        * no arguments: create an empty fixture
        * one argument: host only is specified
        * two arguments: paths to both connection and dependency configuration files are given
        In the latter case, and only in the latter case, unit test data are retrieved.
        """
        
        # Configure complete flag
        conf = connection!=None and dependencies!=None
        if conf:
            assert isinstance(connection, str), "'connection' should be a string"
            assert isinstance(connection, str), "'dependencies' should be a string"
        
        # Cannot specify both host and config file
        assert (host==None or (connection==None and dependencies==None)), "Incompatible set of arguments"
        
        # Default values
        self._fileMap = {}
        self._host = DataHost.UNKNOWN
        self._overwrite = False
        self._distant = ""
        self._local = ""
        
        # Initialize
        if host:
#             assert isinstance(host, DataHost) #TODO only works with Python 3
            self._host = host
        if connection:
            self.configureConnection(connection)
        if conf:
            self.configureDependencies(dependencies)
            self.downloadDependencies()
    
    @staticmethod
    def parseDataHost(name):
        """Parse the name of a data hosting solution (case-insensitive).
        """
        
        assert isinstance(name, str), "name should be a string."
        
        for host in DataHost.__all__: #TODO use iterator
            if host.lower() == name.lower():
                return host
        
        return DataHost.UNKNOWN
    
    def configureConnection(self, connection):
        """Configure connection to host according to a configuration file.
        """
        
        assert isinstance(connection, str), "'connection' should be a string"
        
        # Declare options
        parser = argparse.ArgumentParser()
        parser.add_argument("--host", type=str, help="Hosting solution (only WebDAV is supported for now)")
        parser.add_argument("--overwrite", type=bool, help="Allow overwriting local files if they already exist")
        parser.add_argument("--distant-workspace", type=str, help="Path to distant repository workspace")
        parser.add_argument("--local-workspace", type=str, help="Path to local repository workspace")
        
        # Get config file path
        configFile = PyDataSyncFixture.getExecutablePath()
        configFile = os.path.join(configFile, connection)
        
        # Read config file
        values = []
        with open(configFile) as f:
            for line in f:
                if not line == '\n':
                    values.append("--" + line.replace('\n', ''))
        args = parser.parse_args(values)
        
        # Configure object
        self._host = PyDataSyncFixture.parseDataHost(args.host)
        self._overwrite = args.overwrite
        self._distant = args.distant_workspace
        self._local = args.local_workspace
        
        return True #TODO error case
    
    def configureDependencies(self, dependencies):
        """Read the list of required test files and their alias from the dependency file
        """
        
        assert isinstance(dependencies, str), "'dependencies' should be a string"

        # Get config file path
        filepath = os.path.join(self.getExecutablePath(), dependencies)
        
        #TODO check file exists
        
        # Parse lines
        i = 0
        with open(filepath) as f:
            for line in f:
                self.parseDependency(line)
                i +=1
        
        return i
    
    def parseDependency(self, line, separator='\t'):
        """Parse a line of the dependency file.
        """

        assert isinstance(line, str), "'line' should be a string"
        assert isinstance(separator, str), "'separator' should be a string"
        
        # Remove line breaks
        line = line.replace('\n', '')
        
        # Find alias
        srcDst = line.split(separator)
        srcPath = os.path.join(self._distant, srcDst[0])
        if len(srcDst)>1:
            dstPath = os.path.join(self._local, srcDst[1])
            hasAlias = True
        else:
            dstPath = os.path.join(self._local, srcDst[0])
            hasAlias = False
        
        # Register paths
        self._fileMap[dstPath] = srcPath
        
        return hasAlias
    
    #####################
    ## SYNCHRONIZATION ##
    #####################
    
    def downloadDependencies(self, execute = True):
        """Download files from host.
        """
        
        # Set the number of dependencies
        count = 0
        
        for dst in self._fileMap.keys():
            
            # Get source from destination paths
            src = self._fileMap[dst]
            
            # Download
            
            self.downloadDependency(src, dst, execute)
            count += 1 #TODO fail case
        
        return count
    
    def downloadDependency(self, src, dst, execute = True):
        """Download a file and give it a local name, or alias.
        """
        
        # Build command
        cmd = ""
        
        # iRods
        if self._host == DataHost.IRODS:
            if self._overwrite:
                cmd += "irsync i:"
            else:
                cmd += "iget "
            cmd += src + " " + dst
        
        # Unknown host
        else:
            raise "Unknown data hosting solution."
        
        # Get command
        if not execute:
            return cmd
        
        # Make directory if needed
        dir, file = os.path.split(dst)
        if not os.access(dir, os.F_OK):
            os.mkdir(dir)
        
        # Run command and get output
        #log.info(cmd)
        if not self._overwrite and os.path.exists(dst):
            output = PyDataSyncFixture.executeCommand(cmd)
            return output
    
    def getSourcePath(self, dst):
        """Get the (distant) source path associated to a (local) destination path.
        """
        return self._fileMap[dst]
    
    def getSourcePaths(self):
        """Get the (local) destination path associated to a (distant) source path.
        """
        return self._fileMap.values()
    
    def getDestinationPaths(self):
        """Get the list of all destination paths.
        """
        return self._fileMap.keys()
    
    ################
    ## ATTRIBUTES ##
    ################
    
    #TODO forbid setting member variables
    
    ###########
    ## TOOLS ##
    ###########
    
    @staticmethod
    def executeCommand(cmd):
        """Execute a command and get its output.
        """
        
        # Run command
        p = subprocess.Popen(cmd.split(), \
                             stdout=subprocess.PIPE, \
                             stderr=subprocess.PIPE)
        
        # Read command output
        out, err = p.communicate()
        
        # Return command output
        return out
    
    @staticmethod
    def getExecutablePath():
        """Get the path to the directory of the current executable.
        """
        dir, file = os.path.split(os.path.abspath(__file__))
        return dir
    