""" class to handle algorithm invocation table transactions """
#!pylint: disable=invalid-name

from CBILDataCommon.Util.postgres_dbi import Database


class AlgorithmInvocation(object):
    ''' transaction management for algorithm invocation '''

    def __init__(self, script=None, parameters=None, commit=False, gusConfigFile=None):
        self._database = None
        self._algorithm_invocation_id = None

        self.__get_db_handle(gusConfigFile)
        if script is not None:
            self.insertAlgorithmInvocation(script, parameters, commit)


    def __get_db_handle(self, gusConfigFile):
        ''' create database connection '''
        self._database = Database(gusConfigFile)
        self._database.connect()


    def insertAlgorithmInvocation(self, script, parameters, commit):
        ''' create the entry for the algorithm invocation '''
        sql = """INSERT INTO AlgorithmInvocation
(script_name, script_parameters, commit_mode) VALUES (%s, %s, %s)
RETURNING algorithm_invocation_id"""
        with self._database.cursor() as cursor:
            cursor.execute(sql, (script, parameters, commit))
            self._algorithm_invocation_id = cursor.fetchone()[0]

        self._database.commit()


    def getAlgorithmInvocationId(self):
        ''' return algorithm invocation id '''
        return self._algorithm_invocation_id


    def close(self):
        ''' close db connection '''
        self._database.close()
