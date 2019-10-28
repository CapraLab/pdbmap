#!/usr/env/python3
"""After inclusion of this module, you may open the database with:
with PDBMapSQLdb as db:
   db.do_things...

and the database will close on exit from the 'with' block.

To execute many queries with single connection, prepend your queries with:
    PDBMapSQLdb.open_global_connection()
When the queries are complete, be sure to:
    PDBMapSQLdb.close_global_connection()
"""

import MySQLdb, MySQLdb.cursors
import datetime
import sys
import traceback

import logging
LOGGER = logging.getLogger(__name__)

import gc

from lib.PDBMapGlobals import PDBMapGlobals
assert(PDBMapGlobals.config)


class PDBMapSQLdb(object):
    """Wrap MySQLdb with every needed method, to standardize errors, and confer database independence
    For data access, client callers will normally use
    with PDBMapSQLdb() as db:  # Opens connection to database via dictionary
       ....
       db.execute()
    to ensure proper release of database connections.  Typically the host/user/password parameters to MySQL
    server are accepted from the defaults in the PDBMapGlobalConfig dictionary"""

    _connections_count = 0  # Diagnostic to help track # of SQL connections, which should never exceed 1 in our design
    last_query_time = None # Diagnostic aid to identify connections held open"""

    _static_access_dictionary = {}
    _global_connection = None
    _static_access_dictionary = PDBMapGlobals.config

    @staticmethod
    def set_access_dictionary(global_config_dict):
        """A dictionry with host, user, db, and passwd keys provides a
           convenient global default for database access"""
        for key in ['dbhost','dbuser','dbname','dbpass']:
            assert global_config_dict[key],"Unable to open SQL until key %s defined in the access/config dictionary"%key
        PDBMapSQLdb._static_access_dictionary = global_config_dict 

    @staticmethod
    def open_global_connection(config_dict = None):
        LOGGER.warning("Opening global SQL connection.  Caller must take care to close_global_connection() later")
        if PDBMapSQLdb._global_connection:
            LOGGER.critical("You already have called PDBMapSQLdb.open_global_connection().  Terminating")
            assert PDBMapSQLdb._global_connection == None,"You are attempting to open a global connection when you already have one open"
        PDBMapSQLdb._global_connection = PDBMapSQLdb(config_dict)

    @staticmethod
    def close_global_connection(config_dict = None):
        LOGGER.info("Closing global SQL connection.")
        assert PDBMapSQLdb._global_connection,"No global SQL connection have been opened"
        assert PDBMapSQLdb._connections_count == 1,"There are somehow multiple open SQL connections"
        PDBMapSQLdb._global_connection.__close_cursor_and_connection__()
        PDBMapSQLdb._global_connection = None

    def __init__(self,config_dict = None):
        if not config_dict:
           config_dict = PDBMapSQLdb._static_access_dictionary
        if not config_dict:
            raise Exception("You must provide a dictionary with dbhost/dbuser/dbname/dbpass defined, either here or via PDBMapSQLdb.set_access_dictionary")
        if PDBMapSQLdb._global_connection:
            self._connection_is_global = True
            self._db_connection = PDBMapSQLdb._global_connection
        else:
            self._connection_is_global = False
            self._db_connection = MySQLdb.connect(host = config_dict['dbhost'],user = config_dict['dbuser'],db = config_dict['dbname'],passwd = config_dict['dbpass'])
            PDBMapSQLdb._connections_count += 1
        self._db_cursor = None
        self._init_time = datetime.datetime.now()

    @property
    def rowcount(self):
        return self._db_connection.rowcount

    def activate_row_cursor(self):
        self._db_cursor = self._db_connection.cursor(MySQLdb.cursors.Cursor)

    def activate_dict_cursor(self):
        self._db_cursor = self._db_connection.cursor(MySQLdb.cursors.DictCursor)

    def execute(self, query, params):
        if self._db_cursor is None:
            self.activate_row_cursor()
        if params:
            LOGGER.info("SQL Query: %s"%query%params)
            return self._db_cursor.execute(query, params)
        else:
            LOGGER.info("SQL Query: %s"%query)
            return self._db_cursor.execute(query)

    def fetchone(self):
        if self._db_cursor:
            return self._db_cursor.fetchone()
        return None

    def fetchall(self):
        if self._db_cursor:
            return self._db_cursor.fetchall()
        return None

    def fetchall_yield(self):
        results = self.fetchall()
        if results:
            for result in results:
                yield result

    def __enter__(self):
        # If the user previously called .open_global_connection, then return that existing one
        if PDBMapSQLdb._global_connection:
            return PDBMapSQLdb._global_connection
        return self

    def __close_cursor_and_connection__(self):
        LOGGER.debug("Closing connection  (__close_cursor_and_connection)")
        if self._db_cursor is not None:
            self._db_cursor.close()
            self._db_cursor = None
        if self._db_connection is not None:
            self._db_connection.close() 
            self._db_connection = None
            PDBMapSQLdb._connections_count -= 1

    def __del__(self):
        if not self._connection_is_global:
            self.__close_cursor_and_connection__()

    def __exit__(self,exc_type, exc_value, tb):
        if not self._connection_is_global:
            if exc_type is not None:
                traceback.print_exception(exc_type, exc_value, tb)
            self.__close_cursor_and_connection__()

    @staticmethod
    def all_closed():
        return PDBMapSQLdb._connections_count == 0
        # return self._db_connection is None and self._db_cursor is None

if __name__ == '__main__':
    LOGGER.exception("PDBMapSQLdb.py is a class definition, and should not be called as mainline")
    sys.exit(1)

