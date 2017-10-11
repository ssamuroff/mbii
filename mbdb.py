import MySQLdb as mdb
import numpy as np
from numpy.core.records import fromarrays

class mbdb(object):

    def __init__(self, sqlserver='localhost', user='flanusse', password='mysqlpass',
                 dbname='mb2_hydro', unix_socket='/home/rmandelb.proj/flanusse/mysql/mysql.sock'):
        try:
            self.db = mdb.connect(sqlserver, user, password, dbname, unix_socket=unix_socket)
        except:
            print("Error when connecting to the database")

    def get_columns(self, table):
        """
        Returns the list of columns in the specificied table
        """
        c = self.db.cursor()
        c.execute('SELECT * FROM %s LIMIT 1;'%table)
        return [i[0] for i in c.description]

    def get(self, table, fields, cond=""):
        """
        Returns the sql query as a nice numpy recarray
        expects the list of fields in the format fields='a,b,c'
        """

        sql = "SELECT %s FROM %s WHERE %s;"%(fields, table, cond)
        print sql
        try:
            # prepare a cursor for the query
            cursor = self.db.cursor()
            cursor.execute(sql)
            print("Fetching %d entries" % cursor.rowcount)
        except:
            if fatal_errors:
                    raise
            else:
                print("Error when runnning the SQL command")
                return

        results = fromarrays(np.array(cursor.fetchall()).squeeze().T,names=fields)
        return results

    def get_sql(self, sql, fields):
	try:
            # prepare a cursor for the query
            cursor = self.db.cursor()
            cursor.execute(sql)
            print("Fetching %d entries" % cursor.rowcount)
        except:
            print("Error when runnning the SQL command")
            return
	results = fromarrays(np.array(cursor.fetchall()).squeeze().T,names=fields)
	return results

