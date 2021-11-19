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
import time
import sys
import traceback

import logging
LOGGER = logging.getLogger(__name__)

import gc

from lib.PDBMapGlobals import PDBMapGlobals


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
    # _static_access_dictionary = PDBMapGlobals.config

    @staticmethod
    def set_access_dictionary(global_config_dict):
        """A dictionary with host, user, db, and passwd keys provides a
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
        PDBMapSQLdb.commit()
        PDBMapSQLdb._global_connection.__close_cursor_and_connection__()
        PDBMapSQLdb._global_connection = None

    def __init__(self,config_dict = None):
        if not config_dict:
           config_dict = PDBMapSQLdb._static_access_dictionary
        if not config_dict:
           config_dict = PDBMapGlobals.config
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

    @property
    def description(self):
        return self._db_cursor.description

    def activate_row_cursor(self):
        self._db_cursor = self._db_connection.cursor(MySQLdb.cursors.Cursor)

    def activate_dict_cursor(self):
        self._db_cursor = self._db_connection.cursor(MySQLdb.cursors.DictCursor)

    def set_session_transaction_isolation_read_committed(self):
        """For massive parallel data loads, SQL row locks can become a (huge) problem.
           Calling this function avoids many MariaDB 1205 exceptions"""
        if self._db_cursor is None:
            self.activate_row_cursor()
        self._db_cursor.execute("SET SESSION TRANSACTION ISOLATION LEVEL READ COMMITTED");

    def execute(self, query: str, params=None) -> int:
        """Execute a SQL query, logging the query, and returning rows affected"""
        if self._db_cursor is None:
            self.activate_row_cursor()
        rows_affected = -1
        # If we are logging a CREATE/DROP/SELECT, put that at level WARN
        # Inserts can make a LOT of log entries - so push them to DEBUG level logging
        log_level = logging.DEBUG if query.startswith('INSERT') else logging.INFO
        if params:
            # Worst case is 2 separate strings
            query_log_entry = ""
            try:
                if isinstance(params,list):
                    query_log_entry = query%tuple(params)
                else:
                    query_log_entry = query%params
            except:
                query_log_entry = "Query=%s\nmismatches params suppled=%s"%(query,str(params))
                log_level = logging.CRITICAL  # << We're going to die anyway win execute
        else:
            query_log_entry = query
  
        LOGGER.log(log_level,"SQL Query: %s"%query_log_entry)

        attempts = 0
        retry = True
        while retry:        
            retry = False
            if attempts > 4:
                msg = "Critical - failed to perform database query after 5 attempts"
                LOGGER.critical(msg)
                raise Exception(msg)
 
            try:
                rows_affected = self._db_cursor.execute(query, params)

            except MySQLdb.OperationalError as err: 
                if err.args[0] == 1205: 
                    # Lock wait timeout
                    LOGGER.exception("Operational Error %s"%str(err))
                    time.sleep(5)
                    retry = True
                else:
                    LOGGER.exception("Unhandled operational Error %s"%str(err))
                    raise

            except MySQLdb.DatabaseError as err: 
                LOGGER.exception("Database Error %s"%str(err))
                raise

            except err:
                LOGGER.exception("Unhandled SQL error from execute: %s"%str(err))
                raise  
            else:
                LOGGER.log(log_level,"Successful SQL Query status return: %d 'rows affected'"%rows_affected)
            
            if retry:
                attempts += 1
                LOGGER.warning("Retrying SQL Query - attempt %d"%(attempts+1))
        
        return rows_affected

    def executemany(self,query, sequence_of_params) -> int:
        """Execute a SQL query, logging the query, and returning rows affected"""
        if self._db_cursor is None:
            self.activate_row_cursor()
        rows_affected = -1
        # If we are logging a CREATE/DROP/SELECT, put that at level WARN
        # Inserts can make a LOT of log entries - so push them to DEBUG level logging
        log_level = logging.DEBUG if query.startswith('INSERT') else logging.INFO

        try:
            query_log_entry = "Executemany (count = %d) query %s"%(len(sequence_of_params),sequence_of_params[0]) + "..."
        except:
            query_log_entry = "Query=%s\nmismatches params suppled=%s"%(query,str(sequence_of_params))
            log_level = logging.CRITICAL  # << We're going to die anyway win execute
        else:
            query_log_entry = query_log_entry[0:75] # Truncate this huge string
  
        LOGGER.log(log_level,"SQL Query: %d params %s"%(len(sequence_of_params),query_log_entry))

        attempts = 0
        retry = True
        while retry:        
            retry = False
            if attempts > 4:
                msg = "Critical - failed to perform database query after 5 attempts"
                LOGGER.critical(msg)
                raise Exception(msg)
 
            try:
                rows_affected = self._db_cursor.executemany(query, sequence_of_params)

            except MySQLdb.OperationalError as err: 
                if err.args[0] == 1205: 
                    # Lock wait timeout
                    LOGGER.exception("Operational Error %s"%str(err))
                    time.sleep(5)
                    retry = True
                else:
                    LOGGER.exception("Unhandled operational Error %s"%str(err))
                    raise

            except MySQLdb.DatabaseError as err: 
                LOGGER.exception("Database Error %s"%str(err))
                raise

            except:
                LOGGER.exception("Unhandled SQL error from MySQLdb.execute()")
                raise  
            else:
                LOGGER.log(log_level,"Successful SQL Query status return: %d 'rows affected'"%rows_affected)
            
            if retry:
                attempts += 1
                LOGGER.warning("Retrying SQL Query - attempt %d"%(attempts+1))
        
        return rows_affected


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
        if self._db_cursor is not None:
            LOGGER.debug("Closing mysql cursor  (__close_cursor_and_connection)")
            self._db_cursor.close()
            self._db_cursor = None
        if self._db_connection is not None:
            LOGGER.debug("Closing mysql connection  (__close_cursor_and_connection)")
            self._db_connection.commit()
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

