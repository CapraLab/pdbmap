import pytest
import logging

from lib import PDBMapGlobals
import gc
LOGGER = logging.getLogger()
def test_PDBMapGlobalsLoadFromDefaults():
    config_dict = PDBMapGlobals.config
    
    required_sql_keys = [
       ('dbhost',"Machine where mariadb is running - mysql -h parameter"),
       ('dbuser',"Username for MariaDb - mysql -u parameter"),
       ('dbname',"Default MariaDB Database to USE"),
       ('dbpass',"MariaDB User password")]

    LOGGER.info(str(config_dict))
    for key_info in required_sql_keys:
        assert key_info[0] in config_dict,"%s Entry missing from global config:\n%s"%(key_info[0],key_info[1])

