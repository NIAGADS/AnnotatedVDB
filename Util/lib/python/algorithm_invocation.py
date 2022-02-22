""" class to handle algorithm invocation table transactions """
#!pylint: disable=invalid-name

from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception
from GenomicsDBData.Util.utils import warning

##
# @package util

class AlgorithmInvocation(object):
    """ transaction management for algorithm invocation """

    def __init__(self, script=None, parameters=None, commit=False, gusConfigFile=None):
        self._database = None
        self._algorithm_invocation_id = None

        self.__get_db_handle(gusConfigFile)
        if script is not None:
            self.insert_algorithm_invocation(script, parameters, commit)


    def __get_db_handle(self, gusConfigFile):
        """ create database connection """
        self._database = Database(gusConfigFile)
        self._database.connect()


    def insert_algorithm_invocation(self, script, parameters, commit):
        """ create the entry for the algorithm invocation """
        sql = """INSERT INTO AnnotatedVDB.AlgorithmInvocation
(script_name, script_parameters, commit_mode) VALUES (%s, %s, %s)
RETURNING algorithm_invocation_id"""
        try:
            with self._database.cursor() as cursor:
                cursor.execute(sql, (script, parameters, commit))
                self._algorithm_invocation_id = cursor.fetchone()[0]

            self._database.commit()
            
        except Exception as err:
            raise_pg_exception(err)
        


    def getAlgorithmInvocationId(self):
        """ return algorithm invocation id """
        return self._algorithm_invocation_id


    def close(self):
        """ close db connection """
        self._database.close()
