import pytest
import logging
LOGGER = logging.getLogger()
from lib import PDBMapGlobals
from lib import PDBMapSQLdb 
import gc

def test_connect():
    LOGGER.info("Test short-term 'with' PDBMapSQLdb() invocation")
    with PDBMapSQLdb() as sql_db:
        assert sql_db is not None,"Fundamental inability to connect to the SQL database"
        sql_db.activate_row_cursor()
        dbNameToUse = PDBMapGlobals.config['dbname']
        sql_db.execute("USE %s;"%dbNameToUse,"")

    sql_db2 = PDBMapSQLdb()
    del sql_db2 

    assert PDBMapSQLdb.all_closed(),"All db connections should be closed, but they are not at this point"

    # Test the global mechanism
    LOGGER.info("Test long term 'global' PDBMapSQLdb() invocation")
    PDBMapSQLdb.open_global_connection()
    with PDBMapSQLdb() as sql_db:
        assert sql_db is not None,"Fundamental inability to connect to the SQL database"
        sql_db.activate_row_cursor()

    assert not PDBMapSQLdb.all_closed(),"The global connection should not yet be closed, but it is closed"
    PDBMapSQLdb.close_global_connection()
    assert PDBMapSQLdb.all_closed(),"The global connection should now be closed - but it is not"
