#!/usr/bin/env python

"""PDBMapGlobals.config is a convenient global handle to a dictionary of filesystem and MariaDB resources
If alternate is provided, the dictionary is initialized from pdbmap_sibling_directory/config/global.config"""

import os
import sys
import logging
import configparser
import pprint
LOGGER = logging.getLogger(__name__)

class _PDBMapGlobals_meta(type):
    """helper class allows class level properties - namely PDBMapGlobals.config"""
    _config = None

    def __init__(cls,*args, **kwargs):
        pass

    @property 
    def config(cls):
        """Return PDBMapGlobals.config dictionary.  Init from default config file if not inited."""
        if not cls._config:
            path_to_this_py_file = os.path.dirname(os.path.abspath(__file__))
            path_to_default_global_config = os.path.join(path_to_this_py_file,"../../config/global.config")
            config = configparser.ConfigParser(allow_no_value=True)
            parsed_file_list = config.read([os.path.abspath(path_to_default_global_config)])
            if len(parsed_file_list) < 1:
                LOGGER.critical("""global.config was not found.  
                                   Either create one in pdbmap/../config/global.config or set PDBMapGlobals.config to dictionary """)
                sys.exit(1)
            cls._config = dict(config.items("Genome_PDB_Mapper"))
            
        return cls._config
        
    @config.setter
    def config(cls, external_dictionary):
        """Set PDBMapGlobals.config to a dictionary at runtime."""
        if cls._config:
            LOGGER.info("Replacing previous PDBMapGlobals definition")
        cls._config = external_dictionary
        LOGGER.info("Setting PDBMapGlobal config from dictionary:\n %s",pprint.pformat(cls._config))

class PDBMapGlobals(metaclass=_PDBMapGlobals_meta):
    """The end-user global class that exposes PDBMapGlobals.config"""
    pass

if __name__ == '__main__':
    LOGGER.exception("PDBMapGlobals.py is a class definition, and should not be called as mainline")
    LOGGER.warning("Nonetheless PDBMapGlobals.config = ",str(PDBMapGlobals.config))
    sys.exit(1)
