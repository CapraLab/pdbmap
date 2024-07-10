#!/usr/bin/env python

"""PDBMapGlobals.config is a convenient global handle to a dictionary of filesystem and MariaDB resources
If alternate is provided, the dictionary is initialized from pdbmap_sibling_directory/config/global.config"""

KEYS_TO_SHROUD = ['dbpass', 'dbhost', 'dbuser']

import os
import sys
import logging
import configparser
import copy

import pprint
LOGGER = logging.getLogger(__name__)

class _PDBMapGlobals_meta(type):
    """helper class allows class level properties - namely PDBMapGlobals.config"""
    _config = None
    _exit_function = sys.exit

    def __init__(cls,*args, **kwargs):
        pass

    @property
    def exit(cls): 
        return cls._exit_function

    @exit.setter
    def exit(cls, exit_function):
        assert callable(exit_function)
        cls._exit_function = exit_function

    @property 
    def config(cls):
        """Return PDBMapGlobals.config dictionary.  Init from default config file if not inited."""
        if not cls._config:
            path_to_this_py_file = os.path.dirname(os.path.abspath(__file__))
            path_to_default_global_config = os.path.join(path_to_this_py_file,"../../config/global.config")
            config = configparser.ConfigParser(allow_no_value=True)
            config_filename_abspath = os.path.abspath(path_to_default_global_config)
            parsed_file_list = config.read([config_filename_abspath])
            if len(parsed_file_list) < 1:
                msg = "global.config was not found in %s."%config_filename_abspath + '\n'
                msg +="Either create one in pdbmap/../config/global.config or set PDBMapGlobals.config to dictionary "
                LOGGER.critical(msg)
                sys.exit(msg)
            cls._config = dict(config.items("Genome_PDB_Mapper"))
            LOGGER.info("%d PDBMapGlobal configuration keys loaded from [Genome_PDB_Mapper] section of %s",len(cls._config),config_filename_abspath)
            
        return cls._config
        
    @config.setter
    def config(cls, external_dictionary):
        """Set PDBMapGlobals.config to a dictionary at runtime."""
        if cls._config:
            LOGGER.info("Replacing previous PDBMapGlobals definition")
        cls._config = external_dictionary
        LOGGER.info("Setting PDBMapGlobal config from dictionary:\n %s",pprint.pformat(PDBMapGlobals.shroud_config_dict(cls._config)))

class PDBMapGlobals(metaclass=_PDBMapGlobals_meta):
    """The end-user global class that exposes PDBMapGlobals.config"""
    @staticmethod
    def shroud_config_dict(_config_dict: dict) -> dict:
        _config_dict_shroud_password = copy.deepcopy(_config_dict)
        for key in KEYS_TO_SHROUD:
            value_to_shroud = _config_dict.get(key, '?')
            _config_dict_shroud_password[key] = '*' * len(value_to_shroud)
        return _config_dict_shroud_password

if __name__ == '__main__':
    LOGGER.exception("PDBMapGlobals.py is a class definition, and should not be called as mainline")
    LOGGER.warning("Nonetheless PDBMapGlobals.config = ",str(PDBMapGlobals.config))
    sys.exit(1)
